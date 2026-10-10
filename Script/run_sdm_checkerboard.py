#!/usr/bin/env python3
"""Run SDM2025 checkerboard tests with optional Gaussian noise.

Run SDM and write results in the template directory. Preserve geometry, units and
settings; disable correction in both stages. See checkerboard_guide_zh.md.
"""

import argparse
import csv
from dataclasses import dataclass
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess

import numpy as np

from make_sdm_checkerboard import read_model


GEOMETRY_COLUMNS = (0, 1, 2, 3, 4, 5, 6, 10, 11)
COORDINATE_TOLERANCE = 5.1e-6  # Native fit files use five decimal degrees.
MODEL_TOLERANCE = 5.1e-5       # Native model files use four decimal places.
COMPLETION_MARKER = "End of computations with SDM2025"
MODEL_METRICS = (
    "vector_rmse_m", "relative_vector_error", "amplitude_correlation",
    "nonzero_amplitude_recovery", "zero_leakage_rms_m", "normalized_leakage",
)


@dataclass
class Dataset:
    filename: str
    settings: str
    weight: float
    fixed_direction: bool


@dataclass
class Template:
    prefix: list
    observation_record: str
    datasets: list
    ns: int
    datunit: float
    nheader: int
    niter: int
    relaxation: str
    smoothing: str
    slipout: str
    gdout: list
    dependencies: list
    original_files: list


def parse_template(path):
    """Read records in sdmgetinp.f order; keep scientific records unchanged."""
    records = [line for line in path.read_text().splitlines()
               if line.strip() and not line.lstrip().startswith("#")]
    fields = [shlex.split(line, comments=False) for line in records]
    iearth, idisc, ns = int(fields[0][0]), int(fields[2][0]), int(fields[3][0])
    if iearth not in (0, 1, 2) or idisc not in (0, 1) or not 1 <= ns <= 99:
        raise ValueError("Invalid Earth/fault settings")

    dependencies, filenames = [], []
    if iearth in (1, 2):
        directory, *names = fields[1]
        suffixes = (".ss", ".ds", ".cl") if iearth == 1 else ("",)
        dependencies.extend(Path(directory) / (name + suffix)
                            for name in names for suffix in suffixes)
        filenames.extend(fields[1])

    position = 4
    for _ in range(ns):
        if idisc == 0:
            dependencies.append(Path(fields[position][4]))
            filenames.append(fields[position][4])
            position += 1
        else:
            position += 3 + int(fields[position + 2][0])
    prefix = records[:position]
    observation_record = records[position]
    tokens = fields[position]
    ngd, datunit, nheader = int(tokens[0]), float(tokens[1].upper().replace("D", "E")), int(tokens[2])
    if ngd <= 0 or not np.isfinite(datunit) or datunit <= 0 or nheader < 0:
        raise ValueError("Invalid observation settings")
    position += 1
    datasets = []
    for _ in range(ngd):
        filename = fields[position][0]
        settings, tokens = records[position + 1], fields[position + 1]
        weight, direction = float(tokens[0].upper().replace("D", "E")), int(tokens[1])
        if not np.isfinite(weight) or weight <= 0 or len(tokens) != (4 if direction == 1 else 2):
            raise ValueError("Invalid dataset weight/direction")
        datasets.append(Dataset(filename, settings, weight, direction == 1))
        position += 2
    tokens = fields[position]
    ndpar = int(tokens[0])
    original_files = [Path(tokens[1])] if ndpar else []
    position += 1 + ndpar  # Skip the correction record and its bounds.
    tokens = fields[position]
    niter, relaxation = int(tokens[0]), tokens[1]
    if niter <= 0:
        raise ValueError("Template niter must be positive")
    dependencies.append(Path(relaxation))
    smoothing = records[position + 1]
    slipout = fields[position + 2][0]
    output_names = [name for row in fields[position + 3:] for name in row]
    if ndpar < 0 or len(output_names) != ngd + bool(ndpar):
        raise ValueError("Invalid correction/output records")
    filenames.extend([relaxation, slipout, *output_names, *(d.filename for d in datasets)])
    # SDM uses character*80 filenames and Fortran quoted strings.
    if any(len(os.fsencode(name)) > 80 or any(c in name for c in "'\"\n\r")
           for name in filenames):
        raise ValueError("Invalid SDM filename")
    original_files.extend(Path(name) for name in output_names)
    return Template(prefix, observation_record, datasets, ns, datunit, nheader,
                    niter, relaxation, smoothing, slipout, output_names[:ngd],
                    dependencies, original_files)


def render_template(template, niter, inputs, slipout, outputs):
    """Write the stage configuration with ndpar=0 and no blank lines."""
    records = template.prefix + [template.observation_record]
    for dataset, filename in zip(template.datasets, inputs):
        records.extend((f"'{filename}'", dataset.settings))
    records.extend(("0", f"{niter} '{template.relaxation}'",
                    template.smoothing, f"'{slipout}'",
                    " ".join(f"'{name}'" for name in outputs)))
    return "\n".join(records) + "\n"


@dataclass
class Observations:
    lines: list
    indices: list
    values: np.ndarray


def read_observations(path, nheader, fixed_direction):
    with path.open(newline="") as stream:
        lines = stream.read().splitlines(keepends=True)
    ncolumns = 4 if fixed_direction else 6
    indices = list(range(nheader, len(lines)))
    if not indices or any(len(lines[i].split()) != ncolumns for i in indices):
        raise ValueError(f"{path}: require header plus nonempty {ncolumns}-column rows")
    values = np.array([lines[i].upper().replace("D", "E").split() for i in indices], dtype=float)
    if not np.all(np.isfinite(values)) or np.any(values[:, 3] <= 0):
        raise ValueError(f"{path}: invalid observation/error")
    return Observations(lines, indices, values)


def read_prediction(path, observations):
    lines = path.read_text().splitlines()
    if len(lines) != len(observations.values) + 1 or any(
            len(line.split()) != 6 for line in lines[1:]):
        raise ValueError(f"{path}: incomplete/malformed prediction output")
    values = np.array([line.upper().replace("D", "E").split() for line in lines[1:]], dtype=float)
    if (not np.all(np.isfinite(values)) or np.any(values[:, 5] != 0)
            or not np.allclose(values[:, :2], observations.values[:, :2],
                               rtol=0, atol=COORDINATE_TOLERANCE)):
        raise ValueError(f"{path}: invalid or misaligned prediction")
    return values[:, 4], [line.split()[4] for line in lines[1:]]


def write_synthetic(path, observations, prediction, tokens, noise_scale, rng):
    """Replace only column three, including preserving the original whitespace."""
    synthetic = prediction.copy()
    if noise_scale > 0:
        synthetic += rng.normal(0.0, noise_scale * observations.values[:, 3])
        tokens = [format(value, ".17g") for value in synthetic]
    if not np.all(np.isfinite(synthetic)):
        raise ValueError("nonfinite synthetic observations/noise")
    lines = observations.lines.copy()
    for index, token in zip(observations.indices, tokens):
        match = list(re.finditer(r"\S+", lines[index]))[2]
        lines[index] = lines[index][:match.start()] + token + lines[index][match.end():]
    with path.open("w", newline="") as stream:
        stream.write("".join(lines))
    return synthetic


def check_geometry(truth, other, description):
    if len(truth) != len(other):
        raise ValueError(f"{description}: segment count mismatch")
    for segment, (expected, actual) in enumerate(zip(truth, other), 1):
        if len(expected) != len(actual) or not np.allclose(
                expected[:, GEOMETRY_COLUMNS], actual[:, GEOMETRY_COLUMNS],
                rtol=0, atol=MODEL_TOLERANCE):
            raise ValueError(f"{description}: segment {segment} geometry/patch order mismatch")


def resolution_metrics(truth, recovered):
    """Area-weighted vector errors and amplitude recovery; undefined -> None."""
    area = truth[:, 5] * truth[:, 6]
    original = np.hypot(truth[:, 7], truth[:, 8])
    restored = np.hypot(recovered[:, 7], recovered[:, 8])
    squared_error = np.sum((recovered[:, 7:9] - truth[:, 7:9]) ** 2, axis=1)
    total = np.sum(area)
    error = np.sum(area * squared_error)
    energy = np.sum(area * original ** 2)
    original_mean = np.sum(area * original) / total
    restored_mean = np.sum(area * restored) / total
    centered_original = original - original_mean
    centered_restored = restored - restored_mean
    variance_original = (0.0 if np.all(original == original[0])
                         else np.sum(area * centered_original ** 2))
    variance_restored = (0.0 if np.all(restored == restored[0])
                         else np.sum(area * centered_restored ** 2))
    correlation = None
    if variance_original > 0 and variance_restored > 0:
        correlation = float(np.clip(
            np.sum(area * centered_original * centered_restored)
            / np.sqrt(variance_original * variance_restored), -1, 1))
    high, low = original > 0, original == 0
    high_sum = np.sum(area[high] * original[high])
    leakage = float(np.sqrt(np.sum(area[low] * restored[low] ** 2)
                            / np.sum(area[low]))) if np.any(low) else None
    return dict(
        patch_count=len(truth), area_km2=float(total),
        vector_rmse_m=float(np.sqrt(error / total)),
        relative_vector_error=float(np.sqrt(error / energy)) if energy > 0 else None,
        amplitude_correlation=correlation,
        nonzero_amplitude_recovery=float(np.sum(area[high] * restored[high])
                                        / high_sum) if high_sum > 0 else None,
        zero_leakage_rms_m=leakage,
        normalized_leakage=(float(leakage / (high_sum / np.sum(area[high])))
                            if leakage is not None and high_sum > 0 else None),
    )


def fit_metrics(datasets, observations, synthetic, predictions, datunit):
    """SDM: w_i^2 = (wfm_j/sum wfm) * sigma_i^-2/sum_j sigma^-2."""
    rows = []
    total_weight = sum(dataset.weight for dataset in datasets)
    joint_mse = 0.0
    for index, (dataset, obs, data, pred) in enumerate(
            zip(datasets, observations, synthetic, predictions), 1):
        # Equivalent inverse-variance weights, avoiding overflow for small sigma.
        weights = (np.min(obs.values[:, 3]) / obs.values[:, 3]) ** 2
        residual_m = (data - pred) * datunit
        mse = float(np.sum(weights * residual_m ** 2) / np.sum(weights))
        joint_mse += dataset.weight / total_weight * mse
        rows.append(dict(scope="dataset", dataset_index=index,
                         dataset_name=dataset.filename, observation_count=len(data),
                         normalized_dataset_weight=dataset.weight / total_weight,
                         displacement_wrms_m=float(np.sqrt(mse))))
    rows.append(dict(scope="all", dataset_index="", dataset_name="ALL",
                     observation_count=sum(len(values) for values in synthetic),
                     normalized_dataset_weight=1.0,
                     displacement_wrms_m=float(np.sqrt(joint_mse))))
    return rows


def check_output_paths(directory, generated, protected):
    """Protect source files before deleting or overwriting test outputs."""
    if len(set(generated)) != len(generated):
        raise ValueError("Duplicate output filenames")
    originals = {path.resolve(strict=True) for path in protected}
    for name in generated:
        path = directory / name
        if (len(os.fsencode(name)) > 80 or any(c in name for c in "'\"\n\r")
                or Path(name).name != name or path.is_symlink() or path.resolve() in originals
                or (path.exists() and (not path.is_file() or path.stat().st_nlink > 1))):
            raise ValueError(f"Unsafe output: {path}")


def run_sdm(executable, directory, config, logfile, outputs):
    """Check the completion banner and outputs: Fortran STOP may return zero."""
    with (directory / logfile).open("w") as stream:
        result = subprocess.run([str(executable)], input=config + "\n", text=True,
                                cwd=directory, stdout=stream, stderr=subprocess.STDOUT)
    output = (directory / logfile).read_text()
    if result.returncode != 0 or COMPLETION_MARKER not in output:
        raise RuntimeError(f"SDM did not complete {config}; inspect {directory / logfile}")
    for name in outputs:
        path = directory / name
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing SDM output: {path}")
    for line in output.splitlines():
        if "Convergence" in line:
            print(f"SDM: {line.strip()}", flush=True)


def run_checkerboard(template_path, truth_path, executable, noise_scale=0.0, seed=0):
    if (not np.isfinite(noise_scale) or noise_scale < 0
            or not isinstance(seed, (int, np.integer)) or seed < 0):
        raise ValueError("require finite noise_scale>=0 and integer seed>=0")
    template_path = template_path.resolve(strict=True)
    truth_path = truth_path.resolve(strict=True)
    executable = executable.resolve(strict=True)
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ValueError(f"SDM executable is not executable: {executable}")
    directory = template_path.parent
    template = parse_template(template_path)
    baseline = (directory / template.slipout).resolve(strict=True)
    dependencies = [(directory / name).resolve(strict=True) for name in template.dependencies]
    # Native reads one factor per physical row and cycles the list as needed.
    relaxation = directory / template.relaxation
    factors = np.array([line.upper().replace("D", "E").split()
                        for line in relaxation.read_text().splitlines()], dtype=float)
    if (factors.ndim != 2 or factors.shape[1] != 1 or factors.size == 0
            or not np.all(np.isfinite(factors) & (factors > 0))):
        raise ValueError("Invalid relaxation factors")
    truth = [values for _, values in read_model(truth_path)[1]]
    if len(truth) != template.ns:
        raise ValueError("Model segment count differs from template")
    baseline_model = [values for _, values in read_model(baseline)[1]]
    check_geometry(truth, baseline_model, "original inversion")
    observation_paths = [(directory / dataset.filename).resolve(strict=True)
                         for dataset in template.datasets]
    observations = [read_observations(path, template.nheader, dataset.fixed_direction)
                    for path, dataset in zip(observation_paths, template.datasets)]
    forward_config = "work_" + template_path.name
    inverse_config = "inv" + forward_config
    forward_model = "work_" + truth_path.name
    inverse_model = "inv" + forward_model
    forward_predictions = ["work_" + Path(name).name for name in template.gdout]
    inverse_predictions = ["inv" + name for name in forward_predictions]
    synthetic_names = ["input_" + Path(dataset.filename).name for dataset in template.datasets]
    model_csv, fit_csv = "invwork_resolution_metrics.csv", "invwork_dataset_fit.csv"
    forward_log, inverse_log = "work_forward.log", "invwork_inverse.log"
    forward_outputs = [forward_model, *forward_predictions, "log_" + forward_model,
                       *(f"s{i:02d}_" + forward_model for i in range(1, template.ns + 1))]
    inverse_outputs = [inverse_model, *inverse_predictions, "log_" + inverse_model,
                       *(f"s{i:02d}_" + inverse_model for i in range(1, template.ns + 1))]
    generated = [forward_config, inverse_config, *forward_outputs, *inverse_outputs,
                 *synthetic_names, forward_log, inverse_log, model_csv, fit_csv]
    protected = [template_path, truth_path, baseline, executable,
                 *dependencies, *observation_paths]
    # Protect the template's original prediction outputs too, if already present.
    protected.extend(directory / name for name in template.original_files
                     if (directory / name).exists())
    check_output_paths(directory, generated, protected)
    # Remove only prechecked test artifacts, so no failed phase can reuse old output.
    for name in generated:
        (directory / name).unlink(missing_ok=True)
    shutil.copyfile(truth_path, directory / forward_model)
    (directory / forward_config).write_text(render_template(
        template, 0, [dataset.filename for dataset in template.datasets],
        forward_model, forward_predictions))
    count = sum(len(obs.values) for obs in observations)
    print(f"{len(observations)} datasets, {count} observations, "
          f"{sum(len(block) for block in truth)} patches; noise_scale={noise_scale:g}, seed={seed}",
          flush=True)
    print(f"Forward: {forward_config}", flush=True)
    run_sdm(executable, directory, forward_config, forward_log, forward_outputs)
    forward_result = [values for _, values in read_model(directory / forward_model)[1]]
    check_geometry(truth, forward_result, "forward model")
    for expected, actual in zip(truth, forward_result):
        if not np.allclose(expected[:, 7:9], actual[:, 7:9], rtol=0, atol=MODEL_TOLERANCE):
            raise ValueError("forward output slip differs from the copied checkerboard truth")
    forward_data = [read_prediction(directory / filename, obs)
                    for filename, obs in zip(forward_predictions, observations)]
    rng = np.random.default_rng(seed)
    synthetic = [write_synthetic(directory / name, obs, values, tokens, noise_scale, rng)
                 for name, obs, (values, tokens) in
                 zip(synthetic_names, observations, forward_data)]
    (directory / inverse_config).write_text(render_template(
        template, template.niter, synthetic_names, inverse_model, inverse_predictions))
    print(f"Inverse: {inverse_config}, niter={template.niter}", flush=True)
    run_sdm(executable, directory, inverse_config, inverse_log, inverse_outputs)
    recovered = [values for _, values in read_model(directory / inverse_model)[1]]
    check_geometry(truth, recovered, "recovered model")
    predictions = [read_prediction(directory / filename, obs)[0]
                   for filename, obs in zip(inverse_predictions, observations)]
    model_rows = [dict(scope="all", segment_index="", **resolution_metrics(
        np.vstack(truth), np.vstack(recovered)))]
    model_rows.extend(dict(scope="segment", segment_index=index,
                           **resolution_metrics(expected, actual))
                      for index, (expected, actual) in enumerate(zip(truth, recovered), 1))
    data_rows = fit_metrics(template.datasets, observations, synthetic, predictions,
                            template.datunit)
    for row in model_rows + data_rows:
        row.update(noise_scale=noise_scale, seed=seed)
    for filename, rows in ((model_csv, model_rows), (fit_csv, data_rows)):
        with (directory / filename).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    for row in model_rows:
        label = "All segments" if row["scope"] == "all" else f"Segment {row['segment_index']}"
        values = "; ".join(f"{key}={row[key]:.6g}" if row[key] is not None
                           else f"{key}=N/A" for key in MODEL_METRICS)
        print(f"{label}: patches={row['patch_count']}, area={row['area_km2']:.6g} km^2; {values}")
    for row in data_rows:
        print(f"{row['dataset_name']}: n={row['observation_count']}, "
              f"displacement WRMS={row['displacement_wrms_m']:.6g} m")
    print(f"Resolution metrics: {directory / model_csv}\nDataset fit: {directory / fit_csv}")
    return model_rows, data_rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("template", type=Path, help="SDM inversion template (.inp)")
    parser.add_argument("checkerboard", type=Path, help="Existing checkerboard truth model")
    parser.add_argument("--sdm", type=Path, default=(Path(__file__).resolve().parents[1]
                        / "SDM2025/SourceCode/sdm2025"), help="SDM2025 executable")
    parser.add_argument("--noise-scale", type=float, default=0.0,
                        help="Gaussian noise multiplier of original sigma (default: 0)")
    parser.add_argument("--seed", type=int, default=0, help="Nonnegative RNG seed (default: 0)")
    args = parser.parse_args()
    run_checkerboard(args.template, args.checkerboard, args.sdm, args.noise_scale, args.seed)
