#!/usr/bin/env python3
import argparse
import csv
import os
import shutil
import tempfile
from multiprocessing import Pool
from pathlib import Path

import h5py

from generate_cache_pdbbind_gridfix import (
    APBS_GRID_POINTS,
    APBS_SPAN_ANGSTROM,
    GRIDFIX_VERSION,
    Featurizer,
    assert_cache_shapes,
    compute_distance_grid,
    compute_hydrophobicity_grid,
    featurize_protein_on_grid,
    make_grid_spec,
    pdb2pqr_wrapper,
    pdb_to_mask,
    prepare_protein_pocket,
    protein_center,
    run_apbs,
    save_atomic_features_separate,
)


PREPARED_ROOT_DEFAULT = "/Users/tevfik/Sandbox/github/PHD/data/external_benchmarks/puresnet_prepared"
OUTPUT_ROOT_DEFAULT = "/Users/tevfik/Sandbox/github/PHD/data/external_benchmarks/puresnet_cache_gridfix_v1"
DATASETS = ("coach420_puresnet", "bu48_puresnet")
MANIFEST_FIELDS = ["dataset", "case", "status", "output_h5", "box_size", "error"]


def parse_args():
    parser = argparse.ArgumentParser(description="Generate gridfix H5 caches for prepared Coach420/BU48 cases.")
    parser.add_argument("--prepared-root", default=PREPARED_ROOT_DEFAULT)
    parser.add_argument("--output-root", default=OUTPUT_ROOT_DEFAULT)
    parser.add_argument("--datasets", nargs="+", choices=DATASETS, default=DATASETS)
    parser.add_argument("--case-list", default=None)
    parser.add_argument("--cases", nargs="*", default=None)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--box-size", type=int, default=36)
    parser.add_argument("--resolution", type=float, default=None)
    parser.add_argument("--nproc", type=int, default=1)
    parser.add_argument("--apbs-bin", default="apbs")
    parser.add_argument("--pdb2pqr-bin", default="pdb2pqr30")
    parser.add_argument("--apbs-timeout", type=int, default=120)
    parser.add_argument("--keep-temp", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--fail-on-error", action="store_true")
    return parser.parse_args()


def read_cases(dataset_root, case_list, explicit_cases, limit):
    if explicit_cases:
        cases = [case.strip() for case in explicit_cases if case.strip()]
    elif case_list:
        with open(case_list) as handle:
            cases = [line.strip() for line in handle if line.strip()]
    else:
        cases = sorted(path.name for path in Path(dataset_root).iterdir() if path.is_dir())
    if limit is not None:
        cases = cases[:limit]
    return cases


def write_manifest_header(path):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS)
        writer.writeheader()


def append_manifest(path, row):
    with open(path, "a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS)
        writer.writerow(row)


def require_file(path, description):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        raise RuntimeError(f"Missing or empty {description}: {path}")


def process_benchmark_case(
    dataset,
    case_dir,
    output_h5,
    output_dir,
    case_name,
    box_size,
    resolution,
    apbs_bin,
    pdb2pqr_bin,
    apbs_timeout,
    keep_temp,
    overwrite,
):
    if os.path.exists(output_h5) and not overwrite:
        print(f"[SKIP] Already exists: {output_h5}", flush=True)
        return

    tmp_root = os.path.join(output_dir, "_tmp")
    os.makedirs(tmp_root, exist_ok=True)
    work_dir = tempfile.mkdtemp(prefix=f"{dataset}_{case_name}_", dir=tmp_root)
    tmp_h5 = output_h5 + ".tmp"

    try:
        src_protein_pdb = os.path.join(case_dir, "protein.pdb")
        src_ligand_pdb = os.path.join(case_dir, "ligand.pdb")
        src_fpocket_pdb = os.path.join(case_dir, "fpocket_selected_pocket_atoms.pdb")
        require_file(src_protein_pdb, "protein PDB")
        require_file(src_ligand_pdb, "ligand PDB")
        require_file(src_fpocket_pdb, "selected fpocket pocket")

        protein_pdb = os.path.join(work_dir, "protein.pdb")
        ligand_pdb = os.path.join(work_dir, "ligand.pdb")
        selected_pdb = os.path.join(work_dir, "selected_protein.pdb")
        protein_pqr = os.path.join(work_dir, "selected_protein.pqr")
        calculated_pocket_pdb = os.path.join(work_dir, "calculated_pocket.pdb")
        shutil.copy2(src_protein_pdb, protein_pdb)
        shutil.copy2(src_ligand_pdb, ligand_pdb)

        selected_pdb, calculated_pocket_pdb = prepare_protein_pocket(
            protein_pdb, ligand_pdb, calculated_pocket_pdb, selected_pdb
        )
        pdb2pqr_wrapper(selected_pdb, protein_pqr, pdb2pqr_bin=pdb2pqr_bin)

        center = protein_center(protein_pdb)
        grid = make_grid_spec(center=center, box_size=box_size, resolution=resolution)
        apbs_grid = run_apbs(
            protein_pdb,
            protein_pqr,
            work_dir,
            case_name,
            target=grid,
            apbs_bin=apbs_bin,
            timeout=apbs_timeout,
        )

        featurizer = Featurizer()
        atomic_features = featurize_protein_on_grid(protein_pdb, grid, featurizer)
        feature_arrays = {
            "ligand": pdb_to_mask(ligand_pdb, grid),
            "electrostatic_grid": apbs_grid,
            "shape": pdb_to_mask(protein_pdb, grid),
            "dist_to_ligand": compute_distance_grid(ligand_pdb, grid),
            "hydrophobicity": compute_hydrophobicity_grid(protein_pdb, grid),
            "dist_to_surface": compute_distance_grid(protein_pdb, grid),
        }
        label_arrays = {
            "binding_site_calculated": pdb_to_mask(calculated_pocket_pdb, grid),
            "binding_site_in_dataset": pdb_to_mask(src_fpocket_pdb, grid),
            "binding_site_fpocket_selected": pdb_to_mask(src_fpocket_pdb, grid),
        }
        assert_cache_shapes(feature_arrays, label_arrays, grid)

        os.makedirs(os.path.dirname(output_h5), exist_ok=True)
        with h5py.File(tmp_h5, "w") as h5f:
            save_atomic_features_separate(h5f, atomic_features, featurizer.FEATURE_NAMES)
            for feature_name, feature_array in feature_arrays.items():
                h5f.create_dataset(f"features/{feature_name}", data=feature_array, compression="gzip")
            for label_name, label_array in label_arrays.items():
                h5f.create_dataset(f"label/{label_name}", data=label_array, compression="gzip")
            h5f.attrs["schema_version"] = GRIDFIX_VERSION
            h5f.attrs["benchmark_dataset"] = dataset
            h5f.attrs["box_size"] = box_size
            h5f.attrs["resolution"] = grid.resolution
            h5f.attrs["center"] = grid.center.tolist()
            h5f.attrs["grid_origin"] = grid.origin.tolist()
            h5f.attrs["grid_max"] = grid.max_point.tolist()
            h5f.attrs["physical_span_angstrom"] = grid.span
            h5f.attrs["apbs_grid_points"] = APBS_GRID_POINTS
            h5f.attrs["apbs_span_angstrom"] = APBS_SPAN_ANGSTROM
            h5f.attrs["features"] = ",".join([f"atomic_{name}" for name in featurizer.FEATURE_NAMES] + list(feature_arrays))
            h5f.attrs["labels"] = list(label_arrays)

        os.replace(tmp_h5, output_h5)
        print(f"[OK] {dataset}/{case_name} -> {output_h5}", flush=True)
    except Exception:
        if os.path.exists(tmp_h5):
            os.remove(tmp_h5)
        raise
    finally:
        if not keep_temp:
            shutil.rmtree(work_dir, ignore_errors=True)
            try:
                os.rmdir(tmp_root)
            except OSError:
                pass


def worker(task):
    fail_on_error, args = task[0], task[1:]
    dataset, output_h5, case_name = args[0], args[2], args[4]
    try:
        process_benchmark_case(*args)
        row = {"dataset": dataset, "case": case_name, "status": "ok", "output_h5": output_h5, "box_size": args[5], "error": ""}
    except Exception as exc:
        row = {"dataset": dataset, "case": case_name, "status": "failed", "output_h5": output_h5, "box_size": args[5], "error": repr(exc)}
        print(f"[ERROR] {dataset}/{case_name}: {exc}", flush=True)
        if fail_on_error:
            raise
    return row


def main():
    args = parse_args()
    output_root = Path(args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    manifest_path = output_root / f"manifest_box{args.box_size}.csv"
    write_manifest_header(manifest_path)

    tasks = []
    for dataset in args.datasets:
        dataset_root = Path(args.prepared_root) / dataset
        cases = read_cases(dataset_root, args.case_list, args.cases, args.limit)
        output_dir = output_root / dataset / f"box{args.box_size}"
        for case_name in cases:
            case_dir = dataset_root / case_name
            output_h5 = output_dir / f"{case_name}.h5"
            tasks.append(
                (
                    args.fail_on_error,
                    dataset,
                    str(case_dir),
                    str(output_h5),
                    str(output_dir),
                    case_name,
                    args.box_size,
                    args.resolution,
                    args.apbs_bin,
                    args.pdb2pqr_bin,
                    args.apbs_timeout,
                    args.keep_temp,
                    args.overwrite,
                )
            )

    print(f"Prepared root: {args.prepared_root}")
    print(f"Output root: {output_root}")
    print(f"Tasks: {len(tasks)}")
    print(f"Manifest: {manifest_path}")

    if args.nproc == 1:
        for task in tasks:
            append_manifest(manifest_path, worker(task))
    else:
        with Pool(args.nproc) as pool:
            for row in pool.imap_unordered(worker, tasks):
                append_manifest(manifest_path, row)

    print("ALL DONE.")


if __name__ == "__main__":
    main()
