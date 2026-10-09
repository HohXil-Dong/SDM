#!/usr/bin/env python3
"""Save one 3-D view and one fault-plane view per SDM segment.

SDM separates segments with blank lines in its combined slip output.
Coordinates and patch dimensions are in km; slip amplitude is in m.

Examples (run from the model directory):
    python plot_sdm_slip.py slip_model.dat
    python plot_sdm_slip.py slip_model.dat -o slip.pdf --title "Maduo"
    python plot_sdm_slip.py slip_model.dat --vmax 7 --cmap rainbow --elev 23 --azim -105

The output filename supplies the stem and format: slip.pdf produces
slip_3d.pdf, slip_s01.pdf, slip_s02.pdf, etc. No title is added by default.
"""

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


REQUIRED = (
    "lat_deg", "lon_deg", "depth_km", "x_local_km", "y_local_km",
    "length_km", "width_km", "slp_am_m", "strike_deg", "dip_deg",
)


def read_model(filename):
    """Read the header and the explicit, blank-line-separated SDM segments."""
    with filename.open() as stream:
        names = stream.readline().split()
        missing = [name for name in REQUIRED if name not in names]
        if missing:
            raise ValueError("missing columns: " + ", ".join(missing))
        blocks = [[]]
        for line in stream:
            if line.strip():
                blocks[-1].append(line)
            elif blocks[-1]:
                blocks.append([])

    segments = []
    for block in blocks:
        if not block:
            continue
        segment = np.atleast_1d(np.genfromtxt(block, names=names, dtype=float))
        if any(np.any(~np.isfinite(segment[name])) for name in REQUIRED):
            raise ValueError("the model contains NaN or infinite values")
        if np.any(segment["length_km"] <= 0) or np.any(segment["width_km"] <= 0):
            raise ValueError("patch lengths and widths must be positive")
        if np.any(segment["slp_am_m"] < 0):
            raise ValueError("slip amplitude must be nonnegative")
        segments.append(segment)
    if not segments:
        raise ValueError("the model has no patches")
    return segments


def local_xy(lon, lat):
    """Spherical local east/north in km, relative to min(lon), mean(lat)."""
    radius_km = 6371.0088
    lon0, lat0 = float(np.min(lon)), float(np.mean(lat))
    east = radius_km * np.cos(np.deg2rad(lat0)) * np.deg2rad(lon - lon0)
    north = radius_km * np.deg2rad(lat - lat0)
    return east, north


def corners(east, north, depth, length, width, strike, dip):
    """Centre +/- L/2 along strike +/- W/2 down dip, with elevation = -depth.

    Strike is clockwise from north; down dip follows the right-hand rule.
    Each curved-fault patch is approximated by its local planar rectangle,
    as in the Chen example. Dip angles above 90 degrees remain supported.
    """
    strike, dip = np.deg2rad([strike, dip])
    along = np.array([np.sin(strike), np.cos(strike), 0.0])
    down_dip = np.array([
        np.cos(dip) * np.cos(strike),
        -np.cos(dip) * np.sin(strike),
        -np.sin(dip),
    ])
    centre = np.array([east, north, -depth])
    u, v = 0.5 * length * along, 0.5 * width * down_dip
    return np.array([centre - u - v, centre + u - v,
                     centre + u + v, centre - u + v])


def plot_3d(model, norm, cmap, title, elev, azim):
    east, north = local_xy(model["lon_deg"], model["lat_deg"])
    patches = np.array([
        corners(e, n, row["depth_km"], row["length_km"], row["width_km"],
                row["strike_deg"], row["dip_deg"])
        for e, n, row in zip(east, north, model, strict=True)
    ])
    fig = plt.figure(figsize=(10, 5.5))
    ax = fig.add_subplot(projection="3d")
    poly = Poly3DCollection(
        patches, array=model["slp_am_m"], norm=norm, cmap=cmap,
        edgecolors=(0, 0, 0, 0.18), linewidths=0.18,
    )
    ax.add_collection3d(poly)

    xyz = patches.reshape(-1, 3)
    ranges = np.ptp(xyz, axis=0)
    pad = 0.02 * ranges[0]
    ax.set_xlim(xyz[:, 0].min() - pad, xyz[:, 0].max() + pad)
    ax.set_ylim(xyz[:, 1].min() - pad, xyz[:, 1].max() + pad)
    zticks = np.arange(-10 * np.ceil(-xyz[:, 2].min() / 10), 1, 10)
    ax.set_zticks(zticks, labels=[f"{abs(z):g}" for z in zticks])
    ax.set_zlim(xyz[:, 2].min() - 1, 1)
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(5))
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(3))
    ax.set_xlabel("East (km)", labelpad=10)
    ax.set_ylabel("North (km)", labelpad=10)
    ax.set_zlabel("Depth (km)", labelpad=6)
    ax.set_title(title)
    ax.view_init(elev=elev, azim=azim)
    ax.set_box_aspect((ranges[0], max(ranges[1], 0.12 * ranges[0]), ranges[2]))
    fig.colorbar(poly, ax=ax, orientation="horizontal", shrink=0.65,
                 pad=0.08, label="Slip (m)")
    return fig


def plot_plane(segment, norm, cmap, title):
    """Use the segment's own along-strike/down-dip coordinates, without padding."""
    patches = np.array([
        [(x - length / 2, y - width / 2), (x + length / 2, y - width / 2),
         (x + length / 2, y + width / 2), (x - length / 2, y + width / 2)]
        for x, y, length, width in zip(
            segment["x_local_km"], segment["y_local_km"],
            segment["length_km"], segment["width_km"], strict=True,
        )
    ])
    xy = patches.reshape(-1, 2)
    xmin, ymin = xy.min(axis=0)
    xmax, ymax = xy.max(axis=0)
    length, width = xmax - xmin, ymax - ymin
    size = 8 / max(length, width)
    fig, ax = plt.subplots(figsize=(length * size + 1, width * size + 1.4),
                           layout="constrained")
    poly = PolyCollection(
        patches, array=segment["slp_am_m"], norm=norm, cmap=cmap,
        edgecolors=(0, 0, 0, 0.15), linewidths=0.15,
    )
    ax.add_collection(poly)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymax, ymin)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Along strike (km)")
    ax.set_ylabel("Down dip (km)")
    ax.set_title(title)
    fig.colorbar(poly, ax=ax, orientation="horizontal", aspect=40,
                 pad=0.08, label="Slip (m)")
    return fig


def save_figure(fig, filename):
    fig.savefig(filename, dpi=300, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print(filename)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("filename", nargs="?", type=Path, default=Path("slip_model.dat"),
                        help="SDM model with blank lines separating segments")
    parser.add_argument("-o", "--output", type=Path, help="output filename stem and format")
    parser.add_argument("--title", default="", help="optional title for each figure")
    parser.add_argument("--vmax", type=float,
                        help="colorbar maximum in metres (default: model maximum)")
    parser.add_argument("--cmap", default="viridis", help="Matplotlib colormap")
    parser.add_argument("--elev", type=float, default=23, help="3-D elevation in degrees")
    parser.add_argument("--azim", type=float, default=-105, help="3-D azimuth in degrees")
    args = parser.parse_args()

    segments = read_model(args.filename)
    model = np.concatenate(segments)
    vmax = args.vmax if args.vmax is not None else model["slp_am_m"].max()
    if not np.isfinite(vmax) or vmax <= 0:
        raise ValueError("--vmax must be finite and positive")
    norm = mpl.colors.Normalize(0, vmax, clip=True)
    cmap = mpl.colormaps[args.cmap]
    output = args.output or args.filename.with_suffix(".png")
    output.parent.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({"font.size": 9, "axes.linewidth": 0.8,
                         "xtick.direction": "out", "ytick.direction": "out"})

    fig = plot_3d(model, norm, cmap, args.title, args.elev, args.azim)
    save_figure(fig, output.with_name(f"{output.stem}_3d{output.suffix}"))
    for number, segment in enumerate(segments, 1):
        fig = plot_plane(segment, norm, cmap, args.title)
        save_figure(fig, output.with_name(f"{output.stem}_s{number:02d}{output.suffix}"))


if __name__ == "__main__":
    main()
