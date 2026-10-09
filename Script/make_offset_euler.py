#!/usr/bin/env python3
"""Generate SDM's NEU-offset and geocentric small-rotation matrix.

Usage: python3 make_offset_euler.py stations.dat offset_euler_grns.dat
Input: latitude longitude [deg], north/east positive; # comments allowed.
Output: one header line, then a (3 * nstations, 6) matrix.
Rows: all north, all east, all up; station order follows the input in each block.
Columns: b_N, b_E, b_U [m], omega_X, omega_Y, omega_Z [rad].
Rotation follows u_ECEF = omega cross r, using a sphere of radius 6371000 m.
Offsets are constants in each station's own NEU components, not ECEF translations.
All three observation components must use the same station list and order.
Requires NumPy. Matrix entries are unweighted metres per parameter unit.
"""

import argparse

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


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("stations", help="Two-column latitude/longitude file [deg]")
    parser.add_argument("output", help="Output correction-matrix file")
    args = parser.parse_args()

    stations = np.loadtxt(args.stations, ndmin=2)
    matrix = offset_euler_matrix(stations)
    np.savetxt(
        args.output,
        matrix,
        fmt="%.12e",
        header="rows: N then E then U; columns: "
        "b_N[m] b_E[m] b_U[m] omega_X[rad] omega_Y[rad] omega_Z[rad]",
    )


if __name__ == "__main__":
    main()
