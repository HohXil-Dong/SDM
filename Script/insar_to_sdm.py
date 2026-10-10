#!/usr/bin/env python3
"""Convert LOS or azimuth observations to SDM format.

Usage: python3 Script/insar_to_sdm.py input.txt output.dat
Input:  longitude latitude displacement c_up c_north c_east sigma
The three coefficients are direction cosines, not angles. Select their
one-based columns with -U/-N/-E (defaults: 4/5/6; a permutation of 4/5/6).
Example for ENU input with two header lines: -E 4 -N 5 -U 6 --skiprows 2
Blank lines are ignored; columns after the seventh are ignored.
Output: latitude longitude displacement sigma incidence_deg azimuth_deg
No header is written: use nheader=0 and orientation_flag=0 in SDM.

SDM projects d = u_up*cos(i) + u_north*sin(i)*cos(a)
                           + u_east*sin(i)*sin(a).
Incidence is measured from up (0..180 deg); azimuth clockwise from north
(0..360 deg). Negative c_up therefore gives incidence > 90 deg.
Only the direction of the rounded coefficient vector is retained.
Point order, displacement sign and displacement/error units are preserved.
Set SDM's datunit to metres per input displacement unit. Requires NumPy.
"""

import argparse
from pathlib import Path

import numpy as np


def sdm_observations(data):
    """Convert a seven-column table, already arranged as UNE, to SDM."""
    if data.ndim != 2 or data.shape[1] != 7 or len(data) == 0:
        raise ValueError("Expected a nonempty seven-column observation table")
    if not np.all(np.isfinite(data)):
        raise ValueError("All input values must be finite")
    if np.any(data[:, 6] <= 0):
        raise ValueError("Measurement errors must be positive")
    up, north, east = data[:, 3:6].T
    norm = np.linalg.norm(data[:, 3:6], axis=1)
    # Allow norm error from direction cosines rounded to 0.001.
    if np.any(np.abs(norm - 1.0) > 0.001):
        raise ValueError("Columns 4-6 must be near-unit UNE direction cosines")
    # atan2 recovers the unit direction without clipping or reversing its sign.
    incidence = np.degrees(np.arctan2(np.hypot(north, east), up))
    azimuth = np.degrees(np.arctan2(east, north)) % 360.0
    return np.column_stack((data[:, 1], data[:, 0], data[:, 2], data[:, 6],
                            incidence, azimuth))


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input", type=Path, help="Observation file; error in column 7")
    parser.add_argument("output", type=Path, help="Six-column SDM observation file")
    parser.add_argument("-U", dest="up", type=int, choices=(4, 5, 6), default=4,
                        help="Up direction-cosine column, counted from 1 (default: 4)")
    parser.add_argument("-N", dest="north", type=int, choices=(4, 5, 6), default=5,
                        help="North direction-cosine column, counted from 1 (default: 5)")
    parser.add_argument("-E", dest="east", type=int, choices=(4, 5, 6), default=6,
                        help="East direction-cosine column, counted from 1 (default: 6)")
    parser.add_argument("--skiprows", type=int, default=0, metavar="N",
                        help="Number of physical header lines to skip (default: 0)")
    args = parser.parse_args()
    if len({args.up, args.north, args.east}) != 3:
        parser.error("-U, -N and -E must select different columns")
    if args.skiprows < 0:
        parser.error("--skiprows must be nonnegative")
    if (args.output.resolve() == args.input.resolve()
            or (args.output.exists() and args.output.samefile(args.input))):
        raise ValueError("The output must differ from the input observation file")

    data = np.loadtxt(args.input, ndmin=2, skiprows=args.skiprows,
                      usecols=(0, 1, 2, args.up - 1, args.north - 1, args.east - 1, 6))
    observations = sdm_observations(data)
    np.savetxt(args.output, observations,
               fmt=("%10.6f", "%11.6f", "%12.6f", "%10.6f", "%11.6f", "%11.6f"))
    print(f"{len(observations)} observations written to {args.output}")


if __name__ == "__main__":
    main()
