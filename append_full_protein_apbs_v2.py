#!/usr/bin/env python3
"""Append full-protein APBS v2 channels to existing H5 cache files.

This script computes a new APBS grid from the full protein source structure and
writes electrostatic_grid_v2_full_protein_* channels into the H5 file. It is
source-based and therefore intentionally slower than augment_h5_v2_features.py.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import tempfile
from multiprocessing import Pool
from pathlib import Path

import h5py
import numpy as np
import prody as pr

from feature_schema_metadata import apply_feature_schema_metadata
from augment_h5_v2_features import (
    V2_APBS_PREFIX,
    add_apbs_representations,
    add_cache_derived_features,
    move_feature_to_auxiliary,
    update_feature_attr,
)
from generate_cache_pdbbind_gridfix import (
    GridSpec,
    normalize_pdb_for_apbs,
    pdb2pqr_wrapper,
    run_apbs,
    select_protein_heavy_atoms,
)


AUXILIARY_FEATURES = ("dist_to_ligand", "ligand")
SUMMARY_FIELDS = [
    "status",
    "h5_path",
    "case",
    "source_type",
    "source_protein",
    "created_count",
    "created_features",
    "moved_to_auxiliary",
    "error",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Append full-protein APBS v2 channels to H5 cache files.")
    parser.add_argument("h5_paths", nargs="*", help="One or more H5 files.")
    parser.add_argument("--input-dir", default=None, help="Directory to search for H5 files.")
    parser.add_argument("--glob", default="*.h5", help="Glob used with --input-dir.")
    parser.add_argument("--source-type", choices=("auto", "pdbbind", "external"), default="auto")
    parser.add_argument(
        "--pdbbind-root",
        default="/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS/datasets/pdbbind/refined-set",
    )
    parser.add_argument(
        "--external-prepared-root",
        default="/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS/datasets/external_benchmarks/puresnet_prepared",
    )
    parser.add_argument("--nproc", type=int, default=1)
    parser.add_argument("--apbs-bin", default="apbs")
    parser.add_argument("--pdb2pqr-bin", default="pdb2pqr30")
    parser.add_argument("--apbs-timeout", type=int, default=120)
    parser.add_argument("--full-signed-scale", type=float, default=150.0)
    parser.add_argument("--overwrite", action="store_true", help="Recompute existing v2 APBS channels.")
    parser.add_argument(
        "--keep-leakage-features",
        action="store_true",
        help="Keep dist_to_ligand and ligand under features/. By default they are moved to auxiliary/.",
    )
    parser.add_argument("--summary-csv", default=None, help="CSV summary path.")
    parser.add_argument("--keep-temp", action="store_true")
    parser.add_argument("--fail-on-error", action="store_true")
    return parser.parse_args()


def discover_h5_paths(args: argparse.Namespace) -> list[Path]:
    paths = [Path(path).expanduser().resolve() for path in args.h5_paths]
    if args.input_dir:
        paths.extend(sorted(Path(args.input_dir).expanduser().resolve().rglob(args.glob)))
    unique = []
    seen = set()
    for path in paths:
        if path in seen:
            continue
        seen.add(path)
        if not path.exists():
            raise FileNotFoundError(path)
        if path.suffix != ".h5":
            raise ValueError(f"Expected .h5 file: {path}")
        unique.append(path)
    if not unique:
        raise SystemExit("No H5 files found.")
    return unique


def case_name_from_h5(path: Path) -> str:
    return path.stem


def external_dataset_from_h5(path: Path) -> str | None:
    for part in path.parts:
        if part in {"coach420_puresnet", "bu48_puresnet"}:
            return part
    parent = path.parent.parent.name
    if parent in {"coach420_puresnet", "bu48_puresnet"}:
        return parent
    return None


def find_source_protein(h5_path: Path, args: argparse.Namespace) -> tuple[str, Path]:
    case = case_name_from_h5(h5_path)
    candidates: list[tuple[str, Path]] = []
    if args.source_type in {"auto", "external"}:
        dataset = external_dataset_from_h5(h5_path)
        if dataset:
            candidates.append(
                (
                    "external",
                    Path(args.external_prepared_root).expanduser().resolve() / dataset / case / "protein.pdb",
                )
            )
    if args.source_type in {"auto", "pdbbind"}:
        candidates.append(
            (
                "pdbbind",
                Path(args.pdbbind_root).expanduser().resolve() / case / f"{case}_protein.pdb",
            )
        )

    for source_type, path in candidates:
        if path.exists() and path.stat().st_size > 0:
            return source_type, path
    searched = ", ".join(str(path) for _, path in candidates) or "no candidates"
    raise FileNotFoundError(f"Source protein not found for {h5_path}. Searched: {searched}")


def grid_from_h5(h5f: h5py.File) -> GridSpec:
    box_size = int(h5f.attrs["box_size"])
    resolution = float(h5f.attrs["resolution"])
    origin = np.asarray(h5f.attrs["grid_origin"], dtype=np.float64)
    center = np.asarray(h5f.attrs.get("center", origin + 0.5 * (box_size - 1) * resolution), dtype=np.float64)
    return GridSpec(center=center, box_size=box_size, resolution=resolution, origin=origin)


def prepare_full_protein(source_pdb: Path, output_pdb: Path) -> Path:
    structure = pr.parsePDB(str(source_pdb))
    if structure is None or structure.numAtoms() == 0:
        raise RuntimeError(f"No atoms found in source protein: {source_pdb}")
    protein = select_protein_heavy_atoms(structure, str(source_pdb))
    pr.writePDB(str(output_pdb), protein)
    normalize_pdb_for_apbs(str(output_pdb))
    return output_pdb


def append_one(task: tuple) -> dict[str, str | int]:
    h5_path, args_dict = task
    args = argparse.Namespace(**args_dict)
    h5_path = Path(h5_path)
    case = case_name_from_h5(h5_path)
    tmp_parent = h5_path.parent / "_tmp_v2_full"
    tmp_parent.mkdir(parents=True, exist_ok=True)
    work_dir = tempfile.mkdtemp(prefix=f"{case}_v2_full_", dir=str(tmp_parent))
    created: list[str] = []
    moved: list[str] = []

    try:
        source_type, source_protein = find_source_protein(h5_path, args)
        with h5py.File(h5_path, "r") as h5f:
            grid = grid_from_h5(h5f)
            raw_path = f"features/{V2_APBS_PREFIX}_raw"
            if raw_path in h5f and not args.overwrite:
                return {
                    "status": "skipped",
                    "h5_path": str(h5_path),
                    "case": case,
                    "source_type": source_type,
                    "source_protein": str(source_protein),
                    "created_count": 0,
                    "created_features": "",
                    "moved_to_auxiliary": "",
                    "error": f"{raw_path}_already_exists",
                }

        full_protein_pdb = Path(work_dir) / "full_protein.pdb"
        full_protein_pqr = Path(work_dir) / "full_protein.pqr"
        prepare_full_protein(source_protein, full_protein_pdb)
        pdb2pqr_wrapper(str(full_protein_pdb), str(full_protein_pqr), pdb2pqr_bin=args.pdb2pqr_bin)
        apbs_grid = run_apbs(
            str(full_protein_pdb),
            str(full_protein_pqr),
            work_dir,
            f"{case}_full_protein",
            target=grid,
            apbs_bin=args.apbs_bin,
            timeout=args.apbs_timeout,
        )

        with h5py.File(h5_path, "a") as h5f:
            add_apbs_representations(
                h5f,
                V2_APBS_PREFIX,
                apbs_grid,
                args.overwrite,
                args.full_signed_scale,
                created,
            )
            if not args.keep_leakage_features:
                for feature_name in AUXILIARY_FEATURES:
                    move_feature_to_auxiliary(h5f, feature_name, args.overwrite, moved, [])
            h5f.attrs["v2_full_protein_apbs_status"] = "computed"
            h5f.attrs["v2_full_protein_apbs_source"] = str(source_protein)
            add_cache_derived_features(h5f, args.overwrite, created, args.full_signed_scale)
            update_feature_attr(h5f)
            apply_feature_schema_metadata(h5f)

        return {
            "status": "ok",
            "h5_path": str(h5_path),
            "case": case,
            "source_type": source_type,
            "source_protein": str(source_protein),
            "created_count": len(created),
            "created_features": ";".join(created),
            "moved_to_auxiliary": ";".join(moved),
            "error": "",
        }
    except Exception as exc:
        if args.fail_on_error:
            raise
        return {
            "status": "failed",
            "h5_path": str(h5_path),
            "case": case,
            "source_type": "",
            "source_protein": "",
            "created_count": 0,
            "created_features": "",
            "moved_to_auxiliary": "",
            "error": repr(exc),
        }
    finally:
        if not args.keep_temp:
            shutil.rmtree(work_dir, ignore_errors=True)
            try:
                (h5_path.parent / "_tmp_v2_full").rmdir()
            except OSError:
                pass


def write_row(handle, row: dict[str, str | int]) -> None:
    writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
    writer.writerow(row)
    handle.flush()


def main() -> None:
    args = parse_args()
    h5_paths = discover_h5_paths(args)
    summary_csv = Path(args.summary_csv).expanduser().resolve() if args.summary_csv else h5_paths[0].parent / "v2_full_protein_apbs_summary.csv"
    summary_csv.parent.mkdir(parents=True, exist_ok=True)
    task_args = vars(args).copy()

    total = len(h5_paths)
    completed = ok_count = skipped_count = failed_count = 0

    with open(summary_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()

        def record(row: dict[str, str | int]) -> None:
            nonlocal completed, ok_count, skipped_count, failed_count
            write_row(handle, row)
            completed += 1
            if row["status"] == "ok":
                ok_count += 1
            elif row["status"] == "skipped":
                skipped_count += 1
            elif row["status"] == "failed":
                failed_count += 1
            print(
                f"[PROGRESS] full-protein APBS v2 {completed}/{total} completed | "
                f"remaining={total - completed} | ok={ok_count} | skipped={skipped_count} | "
                f"failed={failed_count} | latest={row['case']} | status={row['status']}",
                flush=True,
            )

        tasks = [(str(path), task_args) for path in h5_paths]
        if args.nproc == 1:
            for task in tasks:
                record(append_one(task))
        else:
            with Pool(processes=args.nproc) as pool:
                for row in pool.imap_unordered(append_one, tasks):
                    record(row)

    print("ALL DONE.", flush=True)
    print(f"OK: {ok_count}", flush=True)
    print(f"Skipped: {skipped_count}", flush=True)
    print(f"Failed: {failed_count}", flush=True)
    print(f"Summary: {summary_csv}", flush=True)
    if failed_count and args.fail_on_error:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
