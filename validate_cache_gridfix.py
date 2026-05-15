import argparse
import csv
import os
from typing import Iterable

import h5py
import numpy as np


CONTINUOUS_FEATURES = {
    "electrostatic_grid",
    "dist_to_ligand",
    "dist_to_surface",
    "hydrophobicity",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Validate generated HDF5 cache tensors and grid metadata.")
    parser.add_argument("--h5-dir", required=True)
    parser.add_argument("--cases", nargs="*", default=None)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--write-csv", default=None)
    parser.add_argument("--fail-on-error", action="store_true")
    return parser.parse_args()


def list_cases(h5_dir: str, explicit_cases: Iterable[str] | None, limit: int | None):
    if explicit_cases:
        cases = [case[:-3] if case.endswith(".h5") else case for case in explicit_cases]
    else:
        cases = sorted(os.path.splitext(name)[0] for name in os.listdir(h5_dir) if name.endswith(".h5"))
    if limit is not None:
        cases = cases[:limit]
    return cases


def validate_one(path: str):
    errors = []
    rows = []
    with h5py.File(path, "r") as h5f:
        box_size = int(h5f.attrs.get("box_size", 0))
        resolution = float(h5f.attrs.get("resolution", np.nan))
        expected_shape = (box_size, box_size, box_size)
        span = resolution * (box_size - 1)
        attr_span = float(h5f.attrs.get("physical_span_angstrom", span))

        if expected_shape[0] <= 0:
            errors.append("missing_or_invalid_box_size")
        if not np.isfinite(resolution) or resolution <= 0:
            errors.append("missing_or_invalid_resolution")
        if not np.isclose(span, attr_span, atol=1e-4):
            errors.append(f"span_attr_mismatch:{span:.6f}!={attr_span:.6f}")

        for group_name in ("features", "label"):
            if group_name not in h5f:
                errors.append(f"missing_group:{group_name}")
                continue
            for name in sorted(h5f[group_name].keys()):
                data = h5f[group_name][name][:]
                if data.shape != expected_shape:
                    errors.append(f"{group_name}/{name}_shape:{data.shape}!={expected_shape}")
                if not np.isfinite(data).all():
                    errors.append(f"{group_name}/{name}_nonfinite")

                nonzero = int(np.count_nonzero(data))
                row = {
                    "file": os.path.basename(path),
                    "group": group_name,
                    "name": name,
                    "shape": "x".join(str(x) for x in data.shape),
                    "nonzero": nonzero,
                    "sum": float(np.sum(data)),
                    "min": float(np.min(data)),
                    "max": float(np.max(data)),
                    "mean": float(np.mean(data)),
                }
                rows.append(row)

                if group_name == "label" and nonzero == 0:
                    errors.append(f"label/{name}_empty")
                if group_name == "features" and name in CONTINUOUS_FEATURES and np.max(data) == np.min(data):
                    errors.append(f"features/{name}_constant")

    return errors, rows


def main():
    args = parse_args()
    cases = list_cases(args.h5_dir, args.cases, args.limit)
    if not cases:
        raise SystemExit(f"No .h5 files found in {args.h5_dir}")

    csv_path = args.write_csv or os.path.join(args.h5_dir, "cache_validation_summary.csv")
    all_rows = []
    failed = {}

    for case in cases:
        path = os.path.join(args.h5_dir, f"{case}.h5")
        if not os.path.exists(path):
            failed[case] = [f"missing_file:{path}"]
            print(f"[FAIL] {case}: missing file")
            continue
        errors, rows = validate_one(path)
        all_rows.extend(rows)
        if errors:
            failed[case] = errors
            print(f"[FAIL] {case}: {'; '.join(errors)}")
        else:
            print(f"[OK] {case}")

    if all_rows:
        with open(csv_path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(all_rows[0].keys()))
            writer.writeheader()
            writer.writerows(all_rows)
        print(f"CSV summary written: {csv_path}")

    print(f"Checked files: {len(cases)}")
    print(f"Failed files: {len(failed)}")
    if failed and args.fail_on_error:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
