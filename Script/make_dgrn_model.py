#!/usr/bin/env python3
"""Extract and plot a local 1D elastic model for DGRN2020.

Usage: python3 Script/make_dgrn_model.py mode lat lon output --model PATH
mode: crust1.0 or ustc2.0; coordinates: latitude, longitude [degrees].
--model: CRUST1 directory containing crust1.bnds/vp/vs/rho, or a USTC
surface-referenced lon/lat/depth/Vp/Vs file. Relative paths use the cwd.
Output: GLOBAL MODEL PARAMETERS block (count, then index/depth/Vp/Vs/rho)
in km, km/s and g/cm^3, plus <output_stem>_velocity.png (Vp red, Vs blue,
positive-down depth). The text block replaces the model section of a
complete DGRN input. Requires NumPy and Matplotlib; see Script/README.md
for the layering, interpolation and empirical density assumptions.
"""

import argparse
from itertools import islice
from pathlib import Path

import numpy as np
from matplotlib.figure import Figure


CRUST1_PARAMETERS = ("bnds", "vp", "vs", "rho")


def normalize_coordinates(lat, lon):
    """Validate degrees and wrap longitude to [-180, 180)."""
    if not np.isfinite(lat) or not np.isfinite(lon) or not -90 <= lat <= 90:
        raise ValueError("Coordinates must be finite; latitude must be in [-90, 90]")
    return lat, (lon + 180) % 360 - 180


def read_crust1(lat, lon, directory):
    """Read the standard CRUST1 files in directory and selected cell centre."""
    lat, lon = normalize_coordinates(lat, lon)
    # Same zero-based indices as getCN1point.f; include the south polar edge.
    ilat = min(int(90 - lat), 179)
    ilon = int(180 + lon)
    row = ilat * 360 + ilon
    values = []
    for parameter in CRUST1_PARAMETERS:
        path = Path(directory) / f"crust1.{parameter}"
        with path.open(encoding="ascii") as source:
            line = next(islice(source, row, row + 1), None)
        if line is None:
            raise ValueError(f"{path}: missing model cell on line {row + 1}")
        record = np.array([float(value) for value in line.split()])
        if record.shape != (9,) or not np.all(np.isfinite(record)):
            raise ValueError(f"{path}: expected nine finite values on line {row + 1}")
        values.append(record)
    return values, (89.5 - ilat, -179.5 + ilon)


def crust1_to_dgrn(boundaries, vp, vs, density):
    """Convert CRUST1 elevations/layer properties to depth/property records."""
    values = np.asarray([boundaries, vp, vs, density], dtype=float)
    if values.shape != (4, 9) or not np.all(np.isfinite(values)):
        raise ValueError("CRUST1 profile must contain four arrays of nine finite values")
    boundaries, vp, vs, density = values
    thickness = boundaries[:-1] - boundaries[1:]
    if np.any(thickness < 0):
        raise ValueError("CRUST1 boundary elevations must be nonincreasing")
    if thickness[0] > 0:
        raise ValueError(
            "This CRUST1 cell has a nonzero-thickness water layer (Vs=0). "
            "DGRN2020 requires a solid elastic model; water cannot be silently removed."
        )
    depths = boundaries[0] - boundaries
    records = []
    for layer in range(8):
        if thickness[layer] == 0:
            continue
        properties = (vp[layer], vs[layer], density[layer])
        top = (depths[layer], *properties)
        # A continuous interface needs only one record.
        if not records or top != records[-1]:
            records.append(top)
        records.append((depths[layer + 1], *properties))
    mantle = (depths[-1], vp[-1], vs[-1], density[-1])
    if not records or mantle != records[-1]:
        records.append(mantle)
    return np.array(records)


def read_ustc(path):
    """Read the observed lon/lat/depth/Vp/Vs grid; require depth/lat/lon order."""
    data = np.loadtxt(path, ndmin=2)
    if data.shape[1] != 5 or not np.all(np.isfinite(data)):
        raise ValueError(f"{path}: expected five finite columns (lon lat depth Vp Vs)")
    lons, lats, depths = (np.unique(data[:, column]) for column in range(3))
    if len(lons) < 2 or len(lats) < 2 or len(depths) < 2:
        raise ValueError(f"{path}: incomplete USTC grid")
    if depths[0] != 0 or np.any(depths < 0):
        raise ValueError(f"{path}: expected nonnegative depths relative to the surface")
    if len(data) != len(lons) * len(lats) * len(depths):
        raise ValueError(f"{path}: missing or extra USTC grid nodes")
    grid = data.reshape(len(depths), len(lats), len(lons), 5)
    if (not np.all(grid[..., 0] == lons[None, None, :])
            or not np.all(grid[..., 1] == lats[None, :, None])
            or not np.all(grid[..., 2] == depths[:, None, None])):
        raise ValueError(f"{path}: grid must be ordered by depth, latitude, longitude")
    return lons, lats, depths, grid[..., 3:5]


def bilinear_profile(lons, lats, velocities, lat, lon):
    """Interpolate Vp/Vs at every depth; velocities have shape (Nz, Ny, Nx, 2)."""
    lat, lon = normalize_coordinates(lat, lon)
    if not lats[0] <= lat <= lats[-1] or not lons[0] <= lon <= lons[-1]:
        raise ValueError(
            f"USTC coordinates outside grid: latitude [{lats[0]:g}, {lats[-1]:g}], "
            f"longitude [{lons[0]:g}, {lons[-1]:g}]; extrapolation is disabled"
        )
    ix = min(int(np.searchsorted(lons, lon, side="right")) - 1, len(lons) - 2)
    iy = min(int(np.searchsorted(lats, lat, side="right")) - 1, len(lats) - 2)
    tx = (lon - lons[ix]) / (lons[ix + 1] - lons[ix])
    ty = (lat - lats[iy]) / (lats[iy + 1] - lats[iy])
    return ((1 - ty) * ((1 - tx) * velocities[:, iy, ix]
                       + tx * velocities[:, iy, ix + 1])
            + ty * ((1 - tx) * velocities[:, iy + 1, ix]
                    + tx * velocities[:, iy + 1, ix + 1]))


def brocher_density(vp):
    """Brocher (2005), eq. (1): Vp [km/s] -> density [g/cm^3], 1.5 <= Vp <= 8.5.

    rho = 1.6612*Vp - 0.4721*Vp^2 + 0.0671*Vp^3
          - 0.0043*Vp^4 + 0.000106*Vp^5.
    Applied to interpolated Vp, including upper-mantle nodes as an explicit
    empirical approximation; it does not supply measured USTC densities.
    """
    vp = np.asarray(vp, dtype=float)
    if not np.all(np.isfinite(vp)) or np.any((vp < 1.5) | (vp > 8.5)):
        raise ValueError("Brocher density requires finite Vp in [1.5, 8.5] km/s")
    return vp * (1.6612 + vp * (-0.4721 + vp * (0.0671 + vp * (-0.0043 + vp * 0.000106))))


def validate_dgrn_model(model):
    """Require a solid layered profile compatible with dgmain/dgsublay."""
    if model.ndim != 2 or model.shape[1] != 4 or len(model) < 2:
        raise ValueError("DGRN2020 requires a layered depth/Vp/Vs/density profile")
    if not np.all(np.isfinite(model)):
        raise ValueError("DGRN model values must be finite")
    depths, vp, vs, density = model.T
    if depths[0] != 0 or np.any(np.diff(depths) < 0) or depths[-1] <= 0:
        raise ValueError("DGRN depths must start at zero and increase towards the half-space")
    if np.any(vp <= 0) or np.any(vs <= 0) or np.any(density <= 0):
        raise ValueError("DGRN2020 requires positive Vp, Vs and density; fluid layers are unsupported")
    if np.any(vp**2 <= 2 * vs**2):
        raise ValueError("DGRN2020 requires lambda > 0, i.e. Vp^2 > 2*Vs^2")
    if np.all(model[:, 1:] == model[0, 1:]):
        raise ValueError("DGRN2020 rejects a homogeneous model (no differential Green's functions)")


def write_dgrn_model(path, model, notes):
    """Write a Fortran-readable model block with # provenance comments."""
    validate_dgrn_model(model)
    lines = ["# DGRN2020 GLOBAL MODEL PARAMETERS", *(f"# {note}" for note in notes),
             "# Depth is positive down from the local surface.",
             "# The last record begins a constant elastic half-space.",
             str(len(model)), "# no depth[km] vp[km/s] vs[km/s] rho[g/cm^3]"]
    for index, row in enumerate(model, 1):
        lines.append(f"{index:4d} " + " ".join(f"{value:.12g}" for value in row))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def velocity_plot_samples(model):
    """Sample the continuous DGRN input profile as depth/Vp/Vs for plotting.

    At each finite interval, rho, rho*Vp^2 (= lambda+2*mu) and rho*Vs^2
    (= mu) vary linearly. Restore speeds from these interpolated quantities.
    Same-depth records preserve the upper/lower sides of a discontinuity.
    Sampling is for drawing only, not DGRN's homogeneous sublayer scheme.
    """
    validate_dgrn_model(model)
    samples = [model[0, :3]]
    for upper, lower in zip(model[:-1], model[1:]):
        if upper[0] == lower[0]:
            samples.append(lower[:3])
            continue
        fraction = np.linspace(0, 1, 101)[1:]
        depths = upper[0] + fraction * (lower[0] - upper[0])
        density = upper[3] + fraction * (lower[3] - upper[3])
        upper_moduli = upper[3] * upper[1:3]**2
        lower_moduli = lower[3] * lower[1:3]**2
        moduli = upper_moduli + fraction[:, None] * (lower_moduli - upper_moduli)
        speeds = np.sqrt(moduli / density[:, None])
        samples.extend(np.column_stack((depths, speeds)))
    return np.array(samples)


def plot_velocity_model(path, model, title):
    """Save the Vp/Vs depth profile without opening a GUI window."""
    samples = velocity_plot_samples(model)
    figure = Figure(figsize=(5.5, 7))
    axes = figure.subplots()
    axes.plot(samples[:, 1], samples[:, 0], color="tab:red", label="Vp", linewidth=1.8)
    axes.plot(samples[:, 2], samples[:, 0], color="tab:blue", label="Vs", linewidth=1.8)
    axes.set_xlabel("Vp & Vs (km/s)")
    axes.set_ylabel("Depth (km, positive down)")
    axes.set_ylim(model[-1, 0], 0)
    axes.set_title(title)
    axes.grid(True, alpha=0.3)
    axes.legend()
    figure.tight_layout()
    figure.savefig(path, dpi=200)
    return figure


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mode", choices=("crust1.0", "ustc2.0"))
    parser.add_argument("lat", type=float, help="latitude [degrees north]")
    parser.add_argument("lon", type=float, help="longitude [degrees east]")
    parser.add_argument("output", type=Path, help="output model block, not a complete DGRN .inp")
    parser.add_argument("--model", type=Path, required=True,
                        help="CRUST1 model directory or USTC surface-referenced model file")
    args = parser.parse_args()
    input_paths = (tuple(args.model / f"crust1.{name}" for name in CRUST1_PARAMETERS)
                   if args.mode == "crust1.0" else (args.model,))
    figure_path = args.output.with_name(args.output.stem + "_velocity.png")
    try:
        lat, lon = normalize_coordinates(args.lat, args.lon)
        protected = {path.resolve() for path in input_paths}
        if any(path.resolve() in protected for path in (args.output, figure_path)):
            raise ValueError("Model and figure outputs must not overwrite a source model file")
        notes = [f"Mode: {args.mode}; requested latitude={lat:.12g}, longitude={lon:.12g} deg",
                 *(f"Source: {path.resolve()}" for path in input_paths)]
        if args.mode == "crust1.0":
            values, (cell_lat, cell_lon) = read_crust1(lat, lon, args.model)
            model = crust1_to_dgrn(*values)
            notes.extend([
                f"Cell centre: latitude={cell_lat:g}, longitude={cell_lon:g} deg; no interpolation",
                f"Surface elevation: {values[0][0]:.12g} km relative to sea level",
                "Original layer properties and mantle retained; zero-thickness layers omitted",
            ])
        else:
            lons, lats, depths, velocities = read_ustc(args.model)
            profile = bilinear_profile(lons, lats, velocities, lat, lon)
            density = brocher_density(profile[:, 0])
            model = np.column_stack((depths, profile, density))
            notes.extend([
                "Vp and Vs interpolated bilinearly in longitude/latitude; depth nodes retained",
                "Density: Brocher (2005) eq. (1), doi:10.1785/0120050077; empirical also in mantle",
                "DGRN uses linear gradients of lambda, mu and density between depth nodes",
            ])
        write_dgrn_model(args.output, model, notes)
        plot_velocity_model(figure_path, model, f"{args.mode}: lat={lat:g}\N{DEGREE SIGN}, lon={lon:g}\N{DEGREE SIGN}")
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(f"Wrote {len(model)} model records to {args.output}; "
          f"half-space begins at {model[-1, 0]:g} km below the local surface.")
    print(f"Saved velocity profile to {figure_path}")


if __name__ == "__main__":
    main()
