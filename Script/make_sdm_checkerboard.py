#!/usr/bin/env python3
r"""Create an SDM checkerboard slip model on the existing fault patches.

Example, with slip in metres and checker cell sizes in kilometres:
    python Script/make_sdm_checkerboard.py slip_model.dat \
        -o checkerboard_9km.dat --cell-size 9 9 \
        --high-slip 1 --low-slip 0 --rake 0

x_local_km is distance along the curved top trace; y_local_km is distance
down dip. Each blank-line-separated segment uses its own local coordinates.
With --origin X Y, a patch centre belongs to the half-open checker cell
[X+i*DX, X+(i+1)*DX) x [Y+j*DY, Y+(j+1)*DY). Even i+j selects high slip;
--phase 1 swaps high and low. The default origin (0, 0) is SDM's native
origin for automatically discretized faults. Boundary cells may be partial;
patches are neither split nor moved to make cell edges fit the patch mesh.

SDM's file convention is slp_strk = S*cos(rake), slp_ddip = -S*sin(rake).
Only these two columns, slp_am_m and rake_deg are replaced. Other tokens,
row order, the header and segment separators are preserved. Stress columns
remain from the input model and MUST be recomputed by SDM forward modeling.
The output must differ from the input. See checkerboard_guide_zh.md beside
this script for the forward and inversion workflow. NumPy is the only dependency.
"""

import argparse
from pathlib import Path
import re

import numpy as np


COLUMNS = (
    "lat_deg", "lon_deg", "depth_km", "x_local_km", "y_local_km",
    "length_km", "width_km", "slp_strk_m", "slp_ddip_m", "slp_am_m",
    "strike_deg", "dip_deg", "rake_deg", "sig_stk_MPa", "sig_ddi_MPa",
    "sig_nrm_MPa", "sig_cmb_MPa",
)
SLIP_COLUMNS = (7, 8, 9, 12)


def read_model(filename):
    """Read the SDM2025 table, keeping source lines and explicit segments.

    Return the original lines and a list of (line_indices, numeric_rows).
    Line indices are zero based; numeric rows remain in their original order.
    The positional schema is checked because SDM reads the first nine
    columns by position, regardless of the column labels.
    """
    lines = filename.read_text().splitlines(keepends=True)
    if not lines or tuple(lines[0].split()) != COLUMNS:
        raise ValueError("expected the 17-column SDM2025 slip-model header")

    blocks = [[]]
    for index, line in enumerate(lines[1:], 1):
        if line.strip():
            if len(line.split()) != len(COLUMNS):
                raise ValueError(f"line {index + 1}: expected 17 columns")
            blocks[-1].append(index)
        elif blocks[-1]:
            blocks.append([])

    segments = []
    for indices in blocks:
        if not indices:
            continue
        model = np.array([lines[index].split() for index in indices], dtype=float)
        if not np.all(np.isfinite(model)):
            raise ValueError("the model contains NaN or infinite values")
        if np.any(model[:, 5:7] <= 0):
            raise ValueError("patch lengths and widths must be positive")
        segments.append((indices, model))
    if not segments:
        raise ValueError("the model has no patches")
    return lines, segments


def make_checkerboard(model, cell_size, high_slip, low_slip, rake,
                      origin=(0.0, 0.0), phase=0):
    """Assign slip from local patch centres, independently of geographic bend.

    Return a modified copy, integer checker indices (strike, dip), and a
    mask selecting high-slip patches. Rake is in degrees; a zero-slip patch
    is labelled rake=0 because its direction is undefined.
    """
    cell_size = np.asarray(cell_size, dtype=float)
    origin = np.asarray(origin, dtype=float)
    if (cell_size.shape != (2,) or not np.all(np.isfinite(cell_size))
            or np.any(cell_size <= 0)):
        raise ValueError("checker cell sizes must be two finite positive values in km")
    if origin.shape != (2,) or not np.all(np.isfinite(origin)):
        raise ValueError("checker origin must be two finite coordinates in km")
    if (not np.all(np.isfinite([high_slip, low_slip]))
            or low_slip < 0 or high_slip <= low_slip):
        raise ValueError("slip amplitudes must satisfy 0 <= low-slip < high-slip")
    if not np.isfinite(rake) or not -180 <= rake <= 360:
        raise ValueError("rake must be finite and between -180 and 360 degrees")
    if phase not in (0, 1):
        raise ValueError("phase must be 0 or 1")

    cells = np.floor((model[:, 3:5] - origin) / cell_size).astype(np.int64)
    high = (cells[:, 0] + cells[:, 1]) % 2 == phase
    amplitude = np.where(high, high_slip, low_slip)
    angle = np.deg2rad(rake)
    updated = model.copy()
    updated[:, 7] = amplitude * np.cos(angle)
    updated[:, 8] = -amplitude * np.sin(angle)
    updated[:, 9] = np.hypot(updated[:, 7], updated[:, 8])
    updated[:, 12] = np.where(amplitude > 0, rake, 0.0)
    return updated, cells, high


def replace_slip_columns(line, row):
    """Replace four numeric tokens; preserve all other text on this line.

    Eight decimals retain sub-micrometre slip precision while keeping the
    whitespace-delimited format accepted by the Fortran free-format reader.
    """
    tokens = list(re.finditer(r"\S+", line))
    for index in reversed(SLIP_COLUMNS):
        token = tokens[index]
        value = row[index]
        replacement = f"{value:.8f}"
        if replacement == "-0.00000000":
            replacement = "0.00000000"
        line = line[:token.start()] + replacement + line[token.end():]
    return line


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("filename", nargs="?", type=Path, default=Path("slip_model.dat"))
    parser.add_argument("-o", "--output", type=Path, default=Path("checkerboard_model.dat"))
    parser.add_argument("--cell-size", type=float, nargs=2, default=(9.0, 9.0),
                        metavar=("STRIKE_KM", "DIP_KM"),
                        help="checker cell sizes along strike and down dip (default: 9 9)")
    parser.add_argument("--high-slip", type=float, default=1.0,
                        help="high-cell slip amplitude in m (default: 1)")
    parser.add_argument("--low-slip", type=float, default=0.0,
                        help="low-cell slip amplitude in m (default: 0)")
    parser.add_argument("--rake", type=float, default=0.0,
                        help="common nonzero-slip rake in degrees (default: 0)")
    parser.add_argument("--origin", type=float, nargs=2, default=(0.0, 0.0),
                        metavar=("X_KM", "Y_KM"),
                        help="checker origin in each segment's local coordinates (default: 0 0)")
    parser.add_argument("--phase", type=int, choices=(0, 1), default=0,
                        help="0: even checker indices have high slip; 1: swap (default: 0)")
    args = parser.parse_args()

    if (args.output.resolve() == args.filename.resolve()
            or (args.output.exists() and args.output.samefile(args.filename))):
        raise ValueError("the output must differ from the input slip model")
    lines, segments = read_model(args.filename)
    reports = []
    for number, (indices, model) in enumerate(segments, 1):
        updated, cells, high = make_checkerboard(
            model, args.cell_size, args.high_slip, args.low_slip, args.rake,
            args.origin, args.phase,
        )
        for index, row in zip(indices, updated):
            lines[index] = replace_slip_columns(lines[index], row)
        nx, ny = (len(np.unique(cells[:, axis])) for axis in (0, 1))
        reports.append(
            f"Segment {number:02d}: {len(model)} patches; {nx} x {ny} checker bins; "
            f"high={np.count_nonzero(high)}, low={np.count_nonzero(~high)}"
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("".join(lines))
    print(f"Saved {args.output}; cells={args.cell_size[0]:g} x {args.cell_size[1]:g} km; "
          f"slip={args.high_slip:g}/{args.low_slip:g} m; rake={args.rake:g} deg")
    for report in reports:
        print(report)
    print("Stress columns are retained from the input; recompute them with SDM forward modeling.")


if __name__ == "__main__":
    main()
