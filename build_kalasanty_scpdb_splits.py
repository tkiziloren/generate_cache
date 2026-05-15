import argparse
import csv
import json
import os
import re
from collections import defaultdict

import h5py
import numpy as np


DEFAULT_H5_DIR = "/Users/tevfik/Sandbox/github/PHD/data/scPDB_cache_gridfix_v1/label_cavity6/box36_span70"
DEFAULT_KALASANTY_DATA_DIR = "data/kalasanty"
DEFAULT_LABEL = "binding_site_cavity6"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build train/validation case lists from generated scPDB H5 caches "
            "using Kalasanty's official folds and exclusions."
        )
    )
    parser.add_argument("--h5-dir", default=DEFAULT_H5_DIR)
    parser.add_argument("--manifest", default=None)
    parser.add_argument("--kalasanty-data-dir", default=DEFAULT_KALASANTY_DATA_DIR)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--label", default=DEFAULT_LABEL)
    parser.add_argument("--folds", type=int, default=10)
    parser.add_argument(
        "--entry-policy",
        choices=["all", "first"],
        default="all",
        help=(
            "Use all scPDB entries for a Kalasanty base PDB id, or only the "
            "first entry. Kalasanty groups entries by base PDB id internally; "
            "our H5 cache stores one file per scPDB entry."
        ),
    )
    return parser.parse_args()


def read_lines(path):
    with open(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def strip_scpdb_suffix(case_id):
    return re.sub(r"_[0-9]+$", "", case_id)


def discover_h5_cases(h5_dir):
    cases = []
    for name in sorted(os.listdir(h5_dir)):
        if name.endswith(".h5"):
            cases.append(name[:-3])
    return cases


def latest_manifest_rows(path):
    if path is None or not os.path.exists(path):
        return {}
    rows = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            case = row.get("case")
            if case:
                rows[case] = row
    return rows


def int_or_none(value):
    if value in (None, ""):
        return None
    return int(float(value))


def label_voxels_from_h5(path, label):
    with h5py.File(path, "r") as h5f:
        if "label" not in h5f or label not in h5f["label"]:
            return None
        return int(np.count_nonzero(h5f["label"][label][:]))


def is_valid_case(case, h5_dir, manifest_rows, label):
    h5_path = os.path.join(h5_dir, f"{case}.h5")
    if not os.path.exists(h5_path):
        return False, "missing_h5"

    row = manifest_rows.get(case, {})
    if row.get("status") == "failed":
        return False, "failed_manifest"

    manifest_label_field = f"{label.replace('binding_site_', '')}_label_voxels"
    label_voxels = int_or_none(row.get(manifest_label_field))
    if label_voxels is None:
        label_voxels = label_voxels_from_h5(h5_path, label)

    if label_voxels is None:
        return False, f"missing_label:{label}"
    if label_voxels <= 0:
        return False, "empty_label"

    atom_field_prefix = label.replace("binding_site_", "")
    atoms = int_or_none(row.get(f"{atom_field_prefix}_label_atoms"))
    atoms_in_box = int_or_none(row.get(f"{atom_field_prefix}_label_atoms_in_box"))
    if atoms is not None and atoms_in_box is not None and atoms_in_box < atoms:
        return False, "label_atoms_outside_box"

    return True, "ok"


def select_cases_for_base_ids(base_ids, cases_by_base, entry_policy):
    selected = []
    for base_id in sorted(base_ids):
        cases = sorted(cases_by_base.get(base_id, []))
        if not cases:
            continue
        if entry_policy == "first":
            selected.append(cases[0])
        else:
            selected.extend(cases)
    return selected


def write_list(path, values):
    with open(path, "w") as handle:
        for value in values:
            handle.write(f"{value}\n")


def main():
    args = parse_args()
    manifest_path = args.manifest or os.path.join(args.h5_dir, "manifest.csv")
    output_dir = args.output_dir or os.path.join(args.h5_dir, "splits_kalasanty")
    os.makedirs(output_dir, exist_ok=True)

    blacklist = set(read_lines(os.path.join(args.kalasanty_data_dir, "scPDB_blacklist.txt")))
    leakage = set(read_lines(os.path.join(args.kalasanty_data_dir, "scPDB_leakage.txt")))
    excluded_cases = blacklist | leakage
    manifest_rows = latest_manifest_rows(manifest_path)
    discovered_cases = discover_h5_cases(args.h5_dir)

    invalid_reasons = defaultdict(int)
    valid_cases = []
    for case in discovered_cases:
        if case in excluded_cases:
            invalid_reasons["kalasanty_excluded"] += 1
            continue
        ok, reason = is_valid_case(case, args.h5_dir, manifest_rows, args.label)
        if ok:
            valid_cases.append(case)
        else:
            invalid_reasons[reason] += 1

    cases_by_base = defaultdict(list)
    for case in valid_cases:
        cases_by_base[strip_scpdb_suffix(case)].append(case)

    valid_base_ids = set(cases_by_base)
    write_list(os.path.join(output_dir, "kalasanty_valid_cases.txt"), sorted(valid_cases))
    write_list(os.path.join(output_dir, "kalasanty_valid_base_ids.txt"), sorted(valid_base_ids))

    fold_summaries = []
    for fold_idx in range(args.folds):
        train_ids = set(read_lines(os.path.join(args.kalasanty_data_dir, f"train_ids_fold{fold_idx}")))
        test_ids = set(read_lines(os.path.join(args.kalasanty_data_dir, f"test_ids_fold{fold_idx}")))

        train_cases = select_cases_for_base_ids(train_ids & valid_base_ids, cases_by_base, args.entry_policy)
        validation_cases = select_cases_for_base_ids(test_ids & valid_base_ids, cases_by_base, args.entry_policy)

        write_list(os.path.join(output_dir, f"fold{fold_idx}_train_cases.txt"), train_cases)
        write_list(os.path.join(output_dir, f"fold{fold_idx}_validation_cases.txt"), validation_cases)
        write_list(os.path.join(output_dir, f"fold{fold_idx}_train_base_ids.txt"), sorted(train_ids & valid_base_ids))
        write_list(os.path.join(output_dir, f"fold{fold_idx}_validation_base_ids.txt"), sorted(test_ids & valid_base_ids))

        fold_summaries.append(
            {
                "fold": fold_idx,
                "train_base_ids": len(train_ids & valid_base_ids),
                "validation_base_ids": len(test_ids & valid_base_ids),
                "train_cases": len(train_cases),
                "validation_cases": len(validation_cases),
                "missing_train_base_ids": len(train_ids - valid_base_ids),
                "missing_validation_base_ids": len(test_ids - valid_base_ids),
            }
        )

    multi_entry_bases = {base: cases for base, cases in cases_by_base.items() if len(cases) > 1}
    with open(os.path.join(output_dir, "summary.json"), "w") as handle:
        json.dump(
            {
                "h5_dir": args.h5_dir,
                "manifest": manifest_path,
                "label": args.label,
                "entry_policy": args.entry_policy,
                "discovered_h5_cases": len(discovered_cases),
                "valid_cases": len(valid_cases),
                "valid_base_ids": len(valid_base_ids),
                "kalasanty_blacklist_cases": len(blacklist),
                "kalasanty_leakage_cases": len(leakage),
                "invalid_reasons": dict(sorted(invalid_reasons.items())),
                "multi_entry_base_ids": len(multi_entry_bases),
                "folds": fold_summaries,
            },
            handle,
            indent=2,
        )

    with open(os.path.join(output_dir, "summary.csv"), "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "fold",
                "train_base_ids",
                "validation_base_ids",
                "train_cases",
                "validation_cases",
                "missing_train_base_ids",
                "missing_validation_base_ids",
            ],
        )
        writer.writeheader()
        writer.writerows(fold_summaries)

    print(f"Output dir: {output_dir}")
    print(f"Discovered H5 cases: {len(discovered_cases)}")
    print(f"Valid cases: {len(valid_cases)}")
    print(f"Valid base PDB ids: {len(valid_base_ids)}")
    print(f"Multi-entry base PDB ids: {len(multi_entry_bases)}")
    print(f"Invalid reasons: {dict(sorted(invalid_reasons.items()))}")


if __name__ == "__main__":
    main()
