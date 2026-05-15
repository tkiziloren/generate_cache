import argparse
import csv
import math
import os
from multiprocessing import Pool

import numpy as np
import prody


DEFAULT_DATASET_ROOT = "/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set"
DEFAULT_SPAN_ANGSTROM = 160.0


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Audit whether PDBBind proteins and pockets fit into the gridfix "
            "point-grid extent. The default box72 gridfix extent is 160 A, "
            "matching the 161-point APBS source grid."
        )
    )
    parser.add_argument("--dataset-root", default=DEFAULT_DATASET_ROOT)
    parser.add_argument("--case-list", default=None)
    parser.add_argument("--box-size", type=int, default=72)
    parser.add_argument("--resolution", type=float, default=None)
    parser.add_argument("--span", type=float, default=DEFAULT_SPAN_ANGSTROM)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--nproc", type=int, default=1)
    parser.add_argument("--write-prefix", default=None)
    return parser.parse_args()


def read_cases(dataset_root, case_list, limit):
    if case_list:
        with open(case_list, "r") as handle:
            cases = [line.strip().split()[0] for line in handle if line.strip()]
    else:
        cases = sorted(
            case
            for case in os.listdir(dataset_root)
            if os.path.isfile(os.path.join(dataset_root, case, f"{case}_protein.pdb"))
        )
    if limit is not None:
        cases = cases[:limit]
    return cases


def coords_from_pdb(path, protein_only=False):
    atoms = prody.parsePDB(path)
    if atoms is None or atoms.numAtoms() == 0:
        return None
    if protein_only:
        protein = atoms.select("protein")
        if protein is not None and protein.numAtoms() > 0:
            atoms = protein
    return atoms.getCoords()


def analyze_case(task):
    case, dataset_root, half_span = task
    case_dir = os.path.join(dataset_root, case)
    protein_pdb = os.path.join(case_dir, f"{case}_protein.pdb")
    pocket_pdb = os.path.join(case_dir, f"{case}_pocket.pdb")
    row = {
        "case": case,
        "status": "ok",
        "full_fit": False,
        "pocket_fit": "",
        "full_required_span": math.nan,
        "full_max_dim": math.nan,
        "pocket_required_span": math.nan,
        "pocket_max_dim": math.nan,
        "error": "",
    }
    try:
        protein_coords = coords_from_pdb(protein_pdb, protein_only=True)
        if protein_coords is None:
            raise RuntimeError(f"No protein atoms in {protein_pdb}")

        center = protein_coords.mean(axis=0)
        protein_max_abs = np.abs(protein_coords - center).max(axis=0)
        protein_dims = protein_coords.max(axis=0) - protein_coords.min(axis=0)
        row["full_required_span"] = float(2.0 * protein_max_abs.max())
        row["full_max_dim"] = float(protein_dims.max())
        row["full_fit"] = bool((protein_max_abs <= half_span + 1e-6).all())

        if os.path.exists(pocket_pdb):
            pocket_coords = coords_from_pdb(pocket_pdb, protein_only=False)
            if pocket_coords is not None:
                pocket_max_abs = np.abs(pocket_coords - center).max(axis=0)
                pocket_dims = pocket_coords.max(axis=0) - pocket_coords.min(axis=0)
                row["pocket_required_span"] = float(2.0 * pocket_max_abs.max())
                row["pocket_max_dim"] = float(pocket_dims.max())
                row["pocket_fit"] = bool((pocket_max_abs <= half_span + 1e-6).all())
        return row
    except Exception as exc:
        row["status"] = "error"
        row["error"] = str(exc)
        return row


def write_outputs(prefix, rows):
    os.makedirs(os.path.dirname(prefix), exist_ok=True) if os.path.dirname(prefix) else None
    csv_path = f"{prefix}_coverage.csv"
    with open(csv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    full_fits = [row["case"] for row in rows if row["status"] == "ok" and row["full_fit"]]
    full_unfits = [
        row for row in rows
        if row["status"] == "ok" and not row["full_fit"]
    ]
    pocket_fits = [
        row["case"] for row in rows
        if row["status"] == "ok" and row["pocket_fit"] is True
    ]
    pocket_unfits = [
        row for row in rows
        if row["status"] == "ok" and row["pocket_fit"] is False
    ]

    with open(f"{prefix}_full_fits.txt", "w") as handle:
        for case in full_fits:
            handle.write(f"{case}\n")
    with open(f"{prefix}_full_unfits.txt", "w") as handle:
        for row in sorted(full_unfits, key=lambda item: item["full_required_span"], reverse=True):
            handle.write(f"{row['case']} requires {row['full_required_span']:.3f}A\n")
    with open(f"{prefix}_pocket_fits.txt", "w") as handle:
        for case in pocket_fits:
            handle.write(f"{case}\n")
    with open(f"{prefix}_pocket_unfits.txt", "w") as handle:
        for row in sorted(pocket_unfits, key=lambda item: item["pocket_required_span"], reverse=True):
            handle.write(f"{row['case']} requires {row['pocket_required_span']:.3f}A\n")

    print(f"Coverage CSV written: {csv_path}")
    print(f"Full-protein fits written: {prefix}_full_fits.txt")
    print(f"Full-protein unfits written: {prefix}_full_unfits.txt")
    print(f"Pocket fits written: {prefix}_pocket_fits.txt")
    print(f"Pocket unfits written: {prefix}_pocket_unfits.txt")


def print_summary(rows, span, box_size, resolution):
    ok_rows = [row for row in rows if row["status"] == "ok"]
    errors = [row for row in rows if row["status"] != "ok"]
    full_fail = [row for row in ok_rows if not row["full_fit"]]
    pocket_known = [row for row in ok_rows if row["pocket_fit"] != ""]
    pocket_fail = [row for row in pocket_known if row["pocket_fit"] is False]

    print(f"Box size: {box_size}")
    print(f"Resolution: {resolution:.8f} A/voxel")
    print(f"Point-grid span: {span:.3f} A")
    print(f"Half span: {span / 2.0:.3f} A")
    print(f"Checked cases: {len(rows)}")
    print(f"Valid cases: {len(ok_rows)}")
    print(f"Errors: {len(errors)}")
    print(f"Full protein fits: {len(ok_rows) - len(full_fail)} / {len(ok_rows)}")
    print(f"Pocket fits: {len(pocket_known) - len(pocket_fail)} / {len(pocket_known)}")

    if full_fail:
        top = sorted(full_fail, key=lambda row: row["full_required_span"], reverse=True)[:15]
        print("Top full-protein failures:")
        for row in top:
            print(f"  {row['case']}: requires {row['full_required_span']:.3f} A")
    if pocket_fail:
        top = sorted(pocket_fail, key=lambda row: row["pocket_required_span"], reverse=True)[:15]
        print("Top pocket failures:")
        for row in top:
            print(f"  {row['case']}: requires {row['pocket_required_span']:.3f} A")


def main():
    args = parse_args()
    prody.confProDy(verbosity="none")
    resolution = args.resolution
    if resolution is None:
        resolution = args.span / float(args.box_size - 1)
    span = resolution * (args.box_size - 1)
    half_span = span / 2.0
    cases = read_cases(args.dataset_root, args.case_list, args.limit)
    tasks = [(case, args.dataset_root, half_span) for case in cases]

    if args.nproc == 1:
        rows = [analyze_case(task) for task in tasks]
    else:
        with Pool(processes=args.nproc) as pool:
            rows = list(pool.imap_unordered(analyze_case, tasks))
        rows.sort(key=lambda row: row["case"])

    print_summary(rows, span, args.box_size, resolution)
    if args.write_prefix:
        write_outputs(args.write_prefix, rows)


if __name__ == "__main__":
    main()
