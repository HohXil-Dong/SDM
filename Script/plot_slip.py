#!/usr/bin/env python3
"""Save slip figures, a patch-corner table and separate segment traces.

SDM separates segments with blank lines in its combined slip output.
Longitude/latitude are in degrees; depth, local coordinates and patch
dimensions are in km; slip amplitude is in m.

Examples (run from the model directory):
    python plot_slip.py slip_model.dat
    python plot_slip.py slip_model.dat -o slip --title "Maduo"
    python plot_slip.py slip_model.dat -o slip --format png
    python plot_slip.py slip_model.dat --vmax 7 --cmap rainbow --elev 23 --azim -105

The output prefix has no image extension: -o slip produces slip_3d.pdf,
slip_s01.pdf, slip_s02.pdf, etc. PDF is the default; --format png selects PNG.
No title is added by default. All segment views share one inches/km scale.
Geometry is saved beside the input model as <model>_patches.dat (12 columns,
UL/UR/LR/LL longitude, latitude, positive-down depth) and <model>_trace.gmt
(two-column longitude/latitude projection of each segment's top patch edges).
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
RADIUS_KM = 6371.0088
PLANE_SIZE_IN = 8.0          # Longest side of the largest segment's drawing area.
COLORBAR_FRACTION = 0.65     # Relative to the drawing area's height.


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
    lon0, lat0 = float(np.min(lon)), float(np.mean(lat))
    east = RADIUS_KM * np.cos(np.deg2rad(lat0)) * np.deg2rad(lon - lon0)
    north = RADIUS_KM * np.deg2rad(lat - lat0)
    return east, north


def corners(east, north, depth, length, width, strike, dip):
    """Centre +/- L/2 along strike +/- W/2 down dip, with elevation = -depth.

    Strike is clockwise from north; down dip follows the right-hand rule.
    Each curved-fault patch is approximated by its local planar rectangle.
    Dip angles above 90 degrees remain supported.
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


def patch_vertices(model):
    """Return UL/UR/LR/LL east, north, elevation corners in model row order."""
    east, north = local_xy(model["lon_deg"], model["lat_deg"])
    return np.array([
        corners(e, n, row["depth_km"], row["length_km"], row["width_km"],
                row["strike_deg"], row["dip_deg"])
        for e, n, row in zip(east, north, model, strict=True)
    ])


def geographic_vertices(model, patches):
    """Invert local_xy: longitude/latitude in degrees, positive-down depth in km."""
    lon0, lat0 = float(np.min(model["lon_deg"])), float(np.mean(model["lat_deg"]))
    vertices = patches.copy()
    vertices[:, :, 0] = lon0 + np.rad2deg(
        patches[:, :, 0] / (RADIUS_KM * np.cos(np.deg2rad(lat0)))
    )
    vertices[:, :, 1] = lat0 + np.rad2deg(patches[:, :, 1] / RADIUS_KM)
    vertices[:, :, 2] = -patches[:, :, 2]
    return vertices


def write_geometry(filename, segments, vertices):
    """Write one 12-column row per patch and independent two-column top traces."""
    patch_file = filename.with_name(f"{filename.stem}_patches.dat")
    trace_file = filename.with_name(f"{filename.stem}_trace.gmt")
    columns = " ".join(
        f"lon_{corner}[deg] lat_{corner}[deg] depth_{corner}[km]"
        for corner in ("UL", "UR", "LR", "LL")
    )
    with patch_file.open("w") as table, trace_file.open("w") as trace:
        table.write(f"# {columns}; depth positive downward\n")
        offset = 0
        for number, segment in enumerate(segments, 1):
            segment_vertices = vertices[offset:offset + len(segment)]
            if number > 1:
                table.write("\n")
            np.savetxt(table, segment_vertices.reshape(-1, 12), fmt="%.10f")

            top = np.flatnonzero(segment["y_local_km"] == segment["y_local_km"].min())
            top = top[np.argsort(segment["x_local_km"][top], kind="stable")]
            trace.write(f"> segment={number}\n")
            np.savetxt(trace, segment_vertices[top, :2, :2].reshape(-1, 2),
                       fmt="%.10f")
            offset += len(segment)
    print(patch_file)
    print(trace_file)


def plot_3d(model, patches, vertices, norm, cmap, title, elev, azim):
    """Draw geographic corners with display proportions set by metric corners."""
    fig = plt.figure(figsize=(10, 5.5))
    ax = fig.add_subplot(projection="3d")
    poly = Poly3DCollection(
        vertices, array=model["slp_am_m"], norm=norm, cmap=cmap,
        edgecolors=(0, 0, 0, 0.18), linewidths=0.18,
    )
    ax.add_collection3d(poly)

    xyz = vertices.reshape(-1, 3)
    # Keep the display proportions in km, even though the horizontal axes use degrees.
    ranges = np.ptp(patches.reshape(-1, 3), axis=0)
    pad_km = 0.02 * max(ranges[:2])
    lat0 = float(np.mean(model["lat_deg"]))
    lon_pad = np.rad2deg(pad_km / (RADIUS_KM * np.cos(np.deg2rad(lat0))))
    lat_pad = np.rad2deg(pad_km / RADIUS_KM)
    ax.set_xlim(xyz[:, 0].min() - lon_pad, xyz[:, 0].max() + lon_pad)
    ax.set_ylim(xyz[:, 1].min() - lat_pad, xyz[:, 1].max() + lat_pad)
    ax.set_zticks(np.arange(0, 10 * np.ceil(xyz[:, 2].max() / 10) + 1, 10))
    ax.set_zlim(xyz[:, 2].max() + 1, -1)
    ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(5))
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(3))
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    ax.set_xlabel("Longitude (°)", labelpad=10)
    ax.set_ylabel("Latitude (°)", labelpad=10)
    ax.set_zlabel("Depth (km)", labelpad=6)
    ax.set_title(title)
    ax.view_init(elev=elev, azim=azim)
    ax.set_box_aspect((max(ranges[0], 0.12 * ranges[1]),
                       max(ranges[1], 0.12 * ranges[0]), max(ranges[2], 2)))
    fig.colorbar(poly, ax=ax, location="right", orientation="vertical",
                 shrink=COLORBAR_FRACTION, pad=0.08, label="Slip (m)")
    return fig


def plane_vertices(segment):
    """Return rectangles in the segment's existing local coordinates."""
    return np.array([
        [(x - length / 2, y - width / 2), (x + length / 2, y - width / 2),
         (x + length / 2, y + width / 2), (x - length / 2, y + width / 2)]
        for x, y, length, width in zip(
            segment["x_local_km"], segment["y_local_km"],
            segment["length_km"], segment["width_km"], strict=True,
        )
    ])


def plot_plane(segment, patches, norm, cmap, title, scale):
    """Draw one segment using the inches/km scale shared by all segments."""
    xy = patches.reshape(-1, 2)
    xmin, ymin = xy.min(axis=0)
    xmax, ymax = xy.max(axis=0)
    length, width = xmax - xmin, ymax - ymin
    plot_width, plot_height = length * scale, width * scale
    # Margins, colorbar gap and width are in inches, independent of segment size.
    left, right, bottom, top = 0.45, 1.0, 0.35, 0.35
    bar_gap, bar_width = 0.2, 0.12
    bar_bottom = bottom + (1 - COLORBAR_FRACTION) / 2 * plot_height
    fig_width = plot_width + (left + right)
    fig_height = plot_height + (bottom + top)
    fig = plt.figure(figsize=(fig_width, fig_height))
    # Fixed physical margins keep text and colorbars from rescaling the patch grid.
    ax = fig.add_axes([left / fig_width, bottom / fig_height,
                       plot_width / fig_width, plot_height / fig_height])
    poly = PolyCollection(
        patches, array=segment["slp_am_m"], norm=norm, cmap=cmap,
        edgecolors=(0, 0, 0, 0.15), linewidths=0.15,
    )
    ax.add_collection(poly)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymax, ymin)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title)
    cax = fig.add_axes([(left + plot_width + bar_gap) / fig_width,
                        bar_bottom / fig_height, bar_width / fig_width,
                        COLORBAR_FRACTION * plot_height / fig_height])
    fig.colorbar(poly, cax=cax, orientation="vertical", label="Slip (m)")
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
    parser.add_argument("-o", "--output", type=Path,
                        help="output filename prefix without an image extension")
    parser.add_argument("--format", choices=("pdf", "png"), default="pdf",
                        help="figure format (default: pdf)")
    parser.add_argument("--title", default="", help="optional title for each figure")
    parser.add_argument("--vmax", type=float,
                        help="colorbar maximum in metres (default: model maximum)")
    parser.add_argument("--cmap", default="viridis", help="Matplotlib colormap")
    parser.add_argument("--elev", type=float, default=23, help="3-D elevation in degrees")
    parser.add_argument("--azim", type=float, default=-105, help="3-D azimuth in degrees")
    args = parser.parse_args()
    if args.output is not None and not args.output.name:
        parser.error("-o expects a nonempty filename prefix")
    if args.output is not None and args.output.suffix.lower() in (".pdf", ".png"):
        parser.error("-o expects a prefix without .pdf/.png; select --format instead")

    segments = read_model(args.filename)
    model = np.concatenate(segments)
    vmax = args.vmax if args.vmax is not None else model["slp_am_m"].max()
    if not np.isfinite(vmax) or vmax <= 0:
        raise ValueError("--vmax must be finite and positive")
    norm = mpl.colors.Normalize(0, vmax, clip=True)
    cmap = mpl.colormaps[args.cmap]
    output = args.output or args.filename.with_suffix("")
    output.parent.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update({"font.size": 9, "axes.linewidth": 0.8,
                         "xtick.direction": "out", "ytick.direction": "out"})

    patches = patch_vertices(model)
    vertices = geographic_vertices(model, patches)
    write_geometry(args.filename, segments, vertices)
    fig = plot_3d(model, patches, vertices, norm, cmap, args.title, args.elev, args.azim)
    save_figure(fig, output.with_name(f"{output.name}_3d.{args.format}"))
    planes = [plane_vertices(segment) for segment in segments]
    scale = PLANE_SIZE_IN / max(np.ptp(plane.reshape(-1, 2), axis=0).max() for plane in planes)
    for number, (segment, plane) in enumerate(zip(segments, planes, strict=True), 1):
        fig = plot_plane(segment, plane, norm, cmap, args.title, scale)
        save_figure(fig, output.with_name(f"{output.name}_s{number:02d}.{args.format}"))


if __name__ == "__main__":
    main()
