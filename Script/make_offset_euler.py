#!/usr/bin/env python3
"""Generate SDM's GPS/InSAR offset and geocentric small-rotation matrix.

Usage: python3 Script/make_offset_euler.py input.dat output.dat --mode gps
       python3 Script/make_offset_euler.py insar_sdm.dat output.dat --mode insar
GPS input: latitude longitude [deg], north/east positive; # comments allowed.
InSAR input: latitude longitude displacement sigma incidence_deg azimuth_deg,
as written by insar_to_sdm.py. Use --inc-col/--azi-col for angle columns
(defaults: 5/6). Incidence is measured from up (0..180 deg); azimuth is
clockwise from north. The projection, including its sign, follows SDM:
G_projected = sin(i)*cos(a)*G_N + sin(i)*sin(a)*G_E + cos(i)*G_U.
Use --lat-col/--lon-col for one-based column numbers (defaults: 1/2),
and --skiprows for header lines (default: 0). Blank lines and other columns
are ignored. The output must differ from the input.
Output: one header line, then a (3 * nstations, 6) GPS matrix or an
(nobservations, 6) InSAR matrix. The default mode is gps.
GPS rows: all north, all east, all up; input order within each block.
InSAR rows: one projected observation per input row, in input order.
Columns: b_N, b_E, b_U [m], omega_X, omega_Y, omega_Z [rad].
Rotation follows u_ECEF = omega cross r, using a sphere of radius 6371000 m.
Offsets are constants in each station's own NEU components, not ECEF translations.
In GPS mode, all three components must use the same station list and order.
Requires NumPy. Matrix entries are unweighted metres per parameter unit.
"""

import argparse
from pathlib import Path

import numpy as np

EARTH_RADIUS_M = 6_371_000.0


def offset_euler_matrix(stations):
    """Build the (3N, 6) matrix from an (N, 2) latitude/longitude array."""
    if stations.ndim != 2 or stations.shape[1] != 2 or len(stations) == 0:
        raise ValueError("Expected a nonempty latitude/longitude table with two columns")
    if not np.all(np.isfinite(stations)):
        raise ValueError("Station coordinates must be finite")
    if np.any(np.abs(stations[:, 0]) > 90):
        raise ValueError("Latitude must be between -90 and 90 degrees")

    lat, lon = np.deg2rad(stations).T
    n = len(stations)
    matrix = np.zeros((3 * n, 6))
    north = matrix[:n]
    east = matrix[n:2 * n]
    up = matrix[2 * n:]

    north[:, 0] = 1.0
    east[:, 1] = 1.0
    up[:, 2] = 1.0

    # Project omega cross r onto each station's north/east/up directions.
    north[:, 3] = EARTH_RADIUS_M * np.sin(lon)
    north[:, 4] = -EARTH_RADIUS_M * np.cos(lon)
    east[:, 3] = -EARTH_RADIUS_M * np.sin(lat) * np.cos(lon)
    east[:, 4] = -EARTH_RADIUS_M * np.sin(lat) * np.sin(lon)
    east[:, 5] = EARTH_RADIUS_M * np.cos(lat)
    # Rotation is tangential to a sphere: the three up-response columns stay zero.
    return matrix


def insar_offset_euler_matrix(stations, angles):
    """Project to each point's (incidence, azimuth) [deg]; return (N, 6)."""
    matrix = offset_euler_matrix(stations)
    n = len(stations)
    if angles.ndim != 2 or angles.shape != (n, 2):
        raise ValueError("Expected one incidence/azimuth pair per observation")
    if not np.all(np.isfinite(angles)):
        raise ValueError("Observation angles must be finite")
    if np.any((angles[:, 0] < 0) | (angles[:, 0] > 180)):
        raise ValueError("Incidence must be between 0 and 180 degrees from up")

    incidence, azimuth = np.deg2rad(angles).T
    c_north = np.sin(incidence) * np.cos(azimuth)
    c_east = np.sin(incidence) * np.sin(azimuth)
    c_up = np.cos(incidence)
    return (c_north[:, None] * matrix[:n]
            + c_east[:, None] * matrix[n:2 * n]
            + c_up[:, None] * matrix[2 * n:])


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("input", type=Path,
                        help="GPS coordinate file or SDM InSAR observation file")
    parser.add_argument("output", type=Path, help="Output correction-matrix file")
    parser.add_argument("--mode", choices=("gps", "insar"), default="gps",
                        help="Observation mode (default: gps)")
    parser.add_argument("--lat-col", type=int, default=1, metavar="COL",
                        help="Latitude column, counted from 1 (default: 1)")
    parser.add_argument("--lon-col", type=int, default=2, metavar="COL",
                        help="Longitude column, counted from 1 (default: 2)")
    parser.add_argument("--inc-col", type=int, default=5, metavar="COL",
                        help="InSAR incidence column [deg], from 1 (default: 5)")
    parser.add_argument("--azi-col", type=int, default=6, metavar="COL",
                        help="InSAR azimuth column [deg], from 1 (default: 6)")
    parser.add_argument("--skiprows", type=int, default=0, metavar="N",
                        help="Number of physical header lines to skip (default: 0)")
    args = parser.parse_args()
    columns = (args.lat_col, args.lon_col)
    if args.mode == "insar":
        columns += (args.inc_col, args.azi_col)
    if min(columns) < 1 or len(set(columns)) != len(columns):
        parser.error("Input columns must be positive and different")
    if args.skiprows < 0:
        parser.error("--skiprows must be nonnegative")
    if (args.output.resolve() == args.input.resolve()
            or (args.output.exists() and args.output.samefile(args.input))):
        raise ValueError("The output must differ from the input coordinate file")

    data = np.loadtxt(args.input, ndmin=2, skiprows=args.skiprows,
                      usecols=tuple(col - 1 for col in columns))
    if args.mode == "gps":
        matrix = offset_euler_matrix(data)
        rows = "N then E then U"
    else:
        matrix = insar_offset_euler_matrix(data[:, :2], data[:, 2:])
        rows = "projected observations in input order"
    np.savetxt(
        args.output,
        matrix,
        fmt="%.12e",
        header=f"rows: {rows}; columns: "
        "b_N[m] b_E[m] b_U[m] omega_X[rad] omega_Y[rad] omega_Z[rad]",
    )


if __name__ == "__main__":
    main()
