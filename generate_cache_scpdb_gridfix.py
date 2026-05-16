import argparse
import csv
import logging
import os
import shutil
import tempfile
from datetime import datetime
from multiprocessing import Pool
from typing import Iterable, Optional

import h5py
import numpy as np
import prody

from feature_schema_metadata import apply_feature_schema_metadata
from generate_cache_pdbbind_gridfix import (
    APBS_GRID_POINTS,
    APBS_SPAN_ANGSTROM,
    GRIDFIX_VERSION,
    Featurizer,
    assert_cache_shapes,
    compute_distance_grid,
    compute_hydrophobicity_grid,
    coords_to_indices,
    featurize_protein_on_grid,
    make_grid_spec,
    mol2_to_pdb,
    pdb2pqr_wrapper,
    pdb_to_mask,
    prepare_protein_pocket,
    protein_center,
    run_apbs,
    save_atomic_features_separate,
)


SCPDB_ROOT_DEFAULT = "/Users/tevfik/Sandbox/github/PHD/data/scPDB"
SCPDB_CACHE_ROOT_DEFAULT = "/Users/tevfik/Sandbox/github/PHD/data/scPDB_cache_gridfix_v1"
DATASET_LABEL_SOURCES = ("site", "cavity6", "cavityALL")
LABEL_DATASET_NAMES = {
    "site": "binding_site_site",
    "cavity6": "binding_site_cavity6",
    "cavityALL": "binding_site_cavityALL",
}
LABEL_MANIFEST_FIELDS = [
    field
    for source in DATASET_LABEL_SOURCES
    for field in (
        f"{source}_label_atoms",
        f"{source}_label_atoms_in_box",
        f"{source}_label_voxels",
    )
]
MANIFEST_FIELDS = [
    "case",
    "status",
    "output_h5",
    "box_size",
    "resolution",
    "target_span",
    "dataset_label_source",
    "protein_atoms",
    "protein_atoms_in_box",
    "ligand_atoms",
    "ligand_atoms_in_box",
    "calculated_label_atoms",
    "calculated_label_atoms_in_box",
    "dataset_label_atoms",
    "dataset_label_atoms_in_box",
    *LABEL_MANIFEST_FIELDS,
    "ligand_voxels",
    "calculated_label_voxels",
    "dataset_label_voxels",
    "shape_voxels",
    "error",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate one HDF5 cache per scPDB case with the same feature/label "
            "schema used by the configurable 3D U-Net training code."
        )
    )
    parser.add_argument("--dataset-root", default=SCPDB_ROOT_DEFAULT)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--case-list", default=None)
    parser.add_argument("--cases", nargs="*", default=None, help="Optional explicit case ids, e.g. 1iki_1 2zb4_1")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--box-size", type=int, default=72)
    parser.add_argument(
        "--target-span",
        type=float,
        default=APBS_SPAN_ANGSTROM,
        help="Physical target grid span in Angstrom. Resolution is target_span / (box_size - 1) unless --resolution is set.",
    )
    parser.add_argument("--resolution", type=float, default=None)
    parser.add_argument(
        "--dataset-label-source",
        choices=DATASET_LABEL_SOURCES,
        default="site",
        help="Source file for label/binding_site_in_dataset.",
    )
    parser.add_argument(
        "--include-all-labels",
        action="store_true",
        help=(
            "Write label/binding_site_site, label/binding_site_cavity6, "
            "and label/binding_site_cavityALL in addition to binding_site_in_dataset."
        ),
    )
    parser.add_argument(
        "--log-file",
        default=None,
        help="Generation log path. Defaults to <output-dir>/generation.log.",
    )
    parser.add_argument("--nproc", type=int, default=1)
    parser.add_argument("--apbs-bin", default="apbs")
    parser.add_argument("--obabel-bin", default="obabel")
    parser.add_argument("--pdb2pqr-bin", default="pdb2pqr30")
    parser.add_argument("--apbs-timeout", type=int, default=120)
    parser.add_argument("--keep-temp", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--fail-on-error", action="store_true")
    return parser.parse_args()


def span_label(target_span: float) -> str:
    if float(target_span).is_integer():
        return str(int(target_span))
    return str(target_span).replace(".", "p")


def default_output_dir(box_size: int, target_span: float, dataset_label_source: str) -> str:
    return os.path.join(
        SCPDB_CACHE_ROOT_DEFAULT,
        f"label_{dataset_label_source}",
        f"box{box_size}_span{span_label(target_span)}",
    )


def read_cases(
    dataset_root: str,
    case_list: Optional[str],
    explicit_cases: Optional[Iterable[str]],
    limit: Optional[int],
):
    if explicit_cases:
        cases = [case.strip() for case in explicit_cases if case.strip()]
    elif case_list:
        with open(case_list, "r") as handle:
            cases = [line.strip() for line in handle if line.strip()]
    else:
        cases = sorted(
            case
            for case in os.listdir(dataset_root)
            if os.path.isdir(os.path.join(dataset_root, case))
        )

    if limit is not None:
        cases = cases[:limit]
    return cases


def scpdb_label_mol2(case_dir: str, dataset_label_source: str) -> str:
    return os.path.join(case_dir, f"{dataset_label_source}.mol2")


def require_file(path: str, description: str):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        raise RuntimeError(f"Missing or empty {description}: {path}")


def required_label_dataset_names(include_all_labels: bool):
    names = ["binding_site_calculated", "binding_site_in_dataset"]
    if include_all_labels:
        names.extend(LABEL_DATASET_NAMES.values())
    return names


def h5_has_required_labels(path: str, include_all_labels: bool):
    try:
        with h5py.File(path, "r") as h5f:
            label_group = h5f.get("label")
            if label_group is None:
                return False
            return all(name in label_group for name in required_label_dataset_names(include_all_labels))
    except OSError:
        return False


def mask_stats(mask: np.ndarray):
    return int(np.count_nonzero(mask)), int(mask.size)


def atom_coverage_stats(pdb_file: str, grid):
    mol = prody.parsePDB(pdb_file)
    if mol is None or mol.numAtoms() == 0:
        return 0, 0
    coords = mol.getCoords()
    _, in_box = coords_to_indices(coords, grid)
    return int(coords.shape[0]), int(np.count_nonzero(in_box))


def write_manifest_header(path: str):
    if os.path.exists(path):
        with open(path, "r", newline="") as handle:
            current_header = next(csv.reader(handle), [])
        if current_header == MANIFEST_FIELDS:
            return
        backup_path = f"{path}.bak_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        os.replace(path, backup_path)
        print(f"[INFO] Existing manifest schema changed. Previous manifest moved to: {backup_path}", flush=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS)
        writer.writeheader()


def append_manifest_row(path: str, row: dict):
    with open(path, "a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS)
        writer.writerow(row)


def blank_label_manifest_values():
    return {field: "" for field in LABEL_MANIFEST_FIELDS}


def setup_logging(log_file: str):
    logger = logging.getLogger("scpdb_gridfix_cache")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    console = logging.StreamHandler()
    console.setFormatter(formatter)
    logger.addHandler(console)

    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    return logger


def log_generation_row(logger, row):
    if row["status"] == "ok":
        logger.info(
            "[OK] %s -> %s | label=%s | dataset_voxels=%s | dataset_atoms_in_box=%s/%s",
            row["case"],
            row["output_h5"],
            row["dataset_label_source"],
            row["dataset_label_voxels"],
            row["dataset_label_atoms_in_box"],
            row["dataset_label_atoms"],
        )
    elif row["status"] == "skipped":
        logger.info("[SKIP] %s -> %s", row["case"], row["output_h5"])
    else:
        logger.error("[FAILED] %s -> %s | %s", row["case"], row["output_h5"], row["error"])


def build_empty_row(
    case_name,
    status,
    output_h5,
    box_size,
    resolution,
    target_span,
    dataset_label_source,
    error="",
):
    return {
        "case": case_name,
        "status": status,
        "output_h5": output_h5,
        "box_size": box_size,
        "resolution": resolution,
        "target_span": target_span,
        "dataset_label_source": dataset_label_source,
        "protein_atoms": "",
        "protein_atoms_in_box": "",
        "ligand_atoms": "",
        "ligand_atoms_in_box": "",
        "calculated_label_atoms": "",
        "calculated_label_atoms_in_box": "",
        "dataset_label_atoms": "",
        "dataset_label_atoms_in_box": "",
        **blank_label_manifest_values(),
        "ligand_voxels": "",
        "calculated_label_voxels": "",
        "dataset_label_voxels": "",
        "shape_voxels": "",
        "error": error,
    }


def process_scpdb_case(
    case_dir,
    output_h5,
    output_dir,
    case_name,
    box_size,
    resolution,
    target_span,
    dataset_label_source,
    apbs_bin,
    obabel_bin,
    pdb2pqr_bin,
    apbs_timeout,
    keep_temp,
    overwrite,
    include_all_labels,
):
    if os.path.exists(output_h5) and not overwrite and h5_has_required_labels(output_h5, include_all_labels):
        return build_empty_row(
            case_name,
            "skipped",
            output_h5,
            box_size,
            resolution,
            target_span,
            dataset_label_source,
        )

    src_protein_mol2 = os.path.join(case_dir, "protein.mol2")
    src_ligand_mol2 = os.path.join(case_dir, "ligand.mol2")
    require_file(src_protein_mol2, "protein mol2")
    require_file(src_ligand_mol2, "ligand mol2")
    label_sources = DATASET_LABEL_SOURCES if include_all_labels else (dataset_label_source,)
    for label_source in label_sources:
        require_file(scpdb_label_mol2(case_dir, label_source), f"{label_source} mol2")

    tmp_root = os.path.join(output_dir, "_tmp")
    os.makedirs(tmp_root, exist_ok=True)
    work_dir = tempfile.mkdtemp(prefix=f"{case_name}_", dir=tmp_root)
    tmp_h5 = output_h5 + ".tmp"

    try:
        protein_pdb = os.path.join(work_dir, "protein.pdb")
        ligand_pdb = os.path.join(work_dir, "ligand.pdb")
        selected_pdb = os.path.join(work_dir, "protein_selected.pdb")
        protein_pqr = os.path.join(work_dir, "protein_selected.pqr")
        calculated_pocket_pdb = os.path.join(work_dir, "binding_site_calculated.pdb")
        label_source_pdbs = {}

        mol2_to_pdb(src_protein_mol2, protein_pdb, obabel_bin=obabel_bin)
        mol2_to_pdb(src_ligand_mol2, ligand_pdb, obabel_bin=obabel_bin)
        for label_source in label_sources:
            label_pdb = os.path.join(work_dir, f"{label_source}.pdb")
            mol2_to_pdb(scpdb_label_mol2(case_dir, label_source), label_pdb, obabel_bin=obabel_bin)
            label_source_pdbs[label_source] = label_pdb

        selected_pdb, calculated_pocket_pdb = prepare_protein_pocket(
            protein_pdb,
            ligand_pdb,
            calculated_pocket_pdb,
            selected_pdb,
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
        dataset_label_pdb = label_source_pdbs[dataset_label_source]
        label_arrays = {
            "binding_site_calculated": pdb_to_mask(calculated_pocket_pdb, grid),
            "binding_site_in_dataset": pdb_to_mask(dataset_label_pdb, grid),
        }
        if include_all_labels:
            for label_source, label_pdb in label_source_pdbs.items():
                label_arrays[LABEL_DATASET_NAMES[label_source]] = pdb_to_mask(label_pdb, grid)

        assert atomic_features.shape[1:] == grid.shape, (
            f"Atomic feature shape {atomic_features.shape[1:]} does not match {grid.shape}"
        )
        assert_cache_shapes(feature_arrays, label_arrays, grid)

        os.makedirs(os.path.dirname(output_h5), exist_ok=True)
        with h5py.File(tmp_h5, "w") as h5f:
            save_atomic_features_separate(h5f, atomic_features, featurizer.FEATURE_NAMES)
            for feature_name, feature_array in feature_arrays.items():
                h5f.create_dataset(f"features/{feature_name}", data=feature_array, compression="gzip")
            for label_name, label_array in label_arrays.items():
                h5f.create_dataset(f"label/{label_name}", data=label_array, compression="gzip")

            h5f.attrs["schema_version"] = f"scpdb_{GRIDFIX_VERSION}"
            h5f.attrs["grid_convention"] = "point_grid_centered_configurable_span"
            h5f.attrs["dataset"] = "scpdb"
            h5f.attrs["case_name"] = case_name
            h5f.attrs["dataset_label_source"] = dataset_label_source
            h5f.attrs["include_all_labels"] = bool(include_all_labels)
            h5f.attrs["all_label_sources"] = list(DATASET_LABEL_SOURCES if include_all_labels else (dataset_label_source,))
            h5f.attrs["box_size"] = box_size
            h5f.attrs["resolution"] = grid.resolution
            h5f.attrs["center"] = grid.center.tolist()
            h5f.attrs["grid_origin"] = grid.origin.tolist()
            h5f.attrs["grid_max"] = grid.max_point.tolist()
            h5f.attrs["physical_span_angstrom"] = grid.span
            h5f.attrs["requested_target_span_angstrom"] = target_span
            h5f.attrs["apbs_grid_points"] = APBS_GRID_POINTS
            h5f.attrs["apbs_span_angstrom"] = APBS_SPAN_ANGSTROM
            h5f.attrs["protein_atom_policy"] = "standard_amino_acid_heavy_atoms_only"
            h5f.attrs["features"] = ",".join(
                [f"atomic_{name}" for name in featurizer.FEATURE_NAMES]
                + list(feature_arrays.keys())
            )
            h5f.attrs["labels"] = list(label_arrays.keys())
            h5f.attrs["electrostatic_grid_shape"] = [box_size, box_size, box_size]
            h5f.attrs["electrostatic_grid_min"] = grid.origin.tolist()
            h5f.attrs["electrostatic_grid_max"] = grid.max_point.tolist()
            apply_feature_schema_metadata(h5f)

        os.replace(tmp_h5, output_h5)

        ligand_voxels, _ = mask_stats(feature_arrays["ligand"])
        calculated_label_voxels, _ = mask_stats(label_arrays["binding_site_calculated"])
        dataset_label_voxels, _ = mask_stats(label_arrays["binding_site_in_dataset"])
        shape_voxels, _ = mask_stats(feature_arrays["shape"])
        protein_atoms, protein_atoms_in_box = atom_coverage_stats(protein_pdb, grid)
        ligand_atoms, ligand_atoms_in_box = atom_coverage_stats(ligand_pdb, grid)
        calculated_label_atoms, calculated_label_atoms_in_box = atom_coverage_stats(calculated_pocket_pdb, grid)
        dataset_label_atoms, dataset_label_atoms_in_box = atom_coverage_stats(dataset_label_pdb, grid)
        label_manifest_values = blank_label_manifest_values()
        for label_source, label_pdb in label_source_pdbs.items():
            source_atoms, source_atoms_in_box = atom_coverage_stats(label_pdb, grid)
            source_dataset_name = LABEL_DATASET_NAMES[label_source]
            if source_dataset_name in label_arrays:
                source_voxels, _ = mask_stats(label_arrays[source_dataset_name])
            else:
                source_voxels, _ = mask_stats(label_arrays["binding_site_in_dataset"])
            label_manifest_values[f"{label_source}_label_atoms"] = source_atoms
            label_manifest_values[f"{label_source}_label_atoms_in_box"] = source_atoms_in_box
            label_manifest_values[f"{label_source}_label_voxels"] = source_voxels

        return {
            "case": case_name,
            "status": "ok",
            "output_h5": output_h5,
            "box_size": box_size,
            "resolution": resolution,
            "target_span": target_span,
            "dataset_label_source": dataset_label_source,
            "protein_atoms": protein_atoms,
            "protein_atoms_in_box": protein_atoms_in_box,
            "ligand_atoms": ligand_atoms,
            "ligand_atoms_in_box": ligand_atoms_in_box,
            "calculated_label_atoms": calculated_label_atoms,
            "calculated_label_atoms_in_box": calculated_label_atoms_in_box,
            "dataset_label_atoms": dataset_label_atoms,
            "dataset_label_atoms_in_box": dataset_label_atoms_in_box,
            **label_manifest_values,
            "ligand_voxels": ligand_voxels,
            "calculated_label_voxels": calculated_label_voxels,
            "dataset_label_voxels": dataset_label_voxels,
            "shape_voxels": shape_voxels,
            "error": "",
        }
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


def process_case(task):
    (
        case_dir,
        output_h5,
        output_dir,
        case_name,
        box_size,
        resolution,
        target_span,
        dataset_label_source,
        apbs_bin,
        obabel_bin,
        pdb2pqr_bin,
        apbs_timeout,
        keep_temp,
        overwrite,
        include_all_labels,
    ) = task
    try:
        return process_scpdb_case(
            case_dir=case_dir,
            output_h5=output_h5,
            output_dir=output_dir,
            case_name=case_name,
            box_size=box_size,
            resolution=resolution,
            target_span=target_span,
            dataset_label_source=dataset_label_source,
            apbs_bin=apbs_bin,
            obabel_bin=obabel_bin,
            pdb2pqr_bin=pdb2pqr_bin,
            apbs_timeout=apbs_timeout,
            keep_temp=keep_temp,
            overwrite=overwrite,
            include_all_labels=include_all_labels,
        )
    except Exception as exc:
        if os.path.exists(output_h5):
            os.remove(output_h5)
        return build_empty_row(
            case_name,
            "failed",
            output_h5,
            box_size,
            resolution,
            target_span,
            dataset_label_source,
            error=str(exc),
        )


def main():
    args = parse_args()
    resolution = args.resolution
    if resolution is None:
        resolution = args.target_span / float(args.box_size - 1)

    output_dir = args.output_dir or default_output_dir(
        args.box_size,
        args.target_span,
        args.dataset_label_source,
    )
    os.makedirs(output_dir, exist_ok=True)

    log_file = args.log_file or os.path.join(output_dir, "generation.log")
    logger = setup_logging(log_file)

    manifest_path = os.path.join(output_dir, "manifest.csv")
    write_manifest_header(manifest_path)

    cases = read_cases(args.dataset_root, args.case_list, args.cases, args.limit)
    logger.info("Gridfix version: scpdb_%s", GRIDFIX_VERSION)
    logger.info("Dataset root: %s", args.dataset_root)
    logger.info("Output dir: %s", output_dir)
    logger.info("Log file: %s", log_file)
    logger.info("Manifest: %s", manifest_path)
    logger.info("Cases: %d", len(cases))
    logger.info("Box size: %d", args.box_size)
    logger.info("Resolution: %.8f A/voxel", resolution)
    logger.info("Target point span: %.3f A", (args.box_size - 1) * resolution)
    logger.info("APBS source: %d points over %.1f A", APBS_GRID_POINTS, APBS_SPAN_ANGSTROM)
    logger.info("Dataset label source: %s.mol2", args.dataset_label_source)
    logger.info("Include all labels: %s", args.include_all_labels)
    if args.include_all_labels:
        logger.info("Extra labels: %s", ", ".join(LABEL_DATASET_NAMES.values()))

    tasks = []
    for case_name in cases:
        case_dir = os.path.join(args.dataset_root, case_name)
        output_h5 = os.path.join(output_dir, f"{case_name}.h5")
        tasks.append(
            (
                case_dir,
                output_h5,
                output_dir,
                case_name,
                args.box_size,
                resolution,
                args.target_span,
                args.dataset_label_source,
                args.apbs_bin,
                args.obabel_bin,
                args.pdb2pqr_bin,
                args.apbs_timeout,
                args.keep_temp,
                args.overwrite,
                args.include_all_labels,
            )
        )

    rows = []
    total_tasks = len(tasks)
    completed = 0
    ok_count = 0
    skipped_count = 0
    failed_count = 0

    def record_progress(row):
        nonlocal completed, ok_count, skipped_count, failed_count
        append_manifest_row(manifest_path, row)
        log_generation_row(logger, row)
        rows.append(row)
        completed += 1
        if row["status"] == "ok":
            ok_count += 1
        elif row["status"] == "skipped":
            skipped_count += 1
        elif row["status"] == "failed":
            failed_count += 1
        remaining = total_tasks - completed
        logger.info(
            "[PROGRESS] scPDB cache %d/%d completed | remaining=%d | ok=%d | skipped=%d | failed=%d | latest=%s",
            completed,
            total_tasks,
            remaining,
            ok_count,
            skipped_count,
            failed_count,
            row["case"],
        )

    if args.nproc == 1:
        for task in tasks:
            row = process_case(task)
            record_progress(row)
    else:
        with Pool(processes=args.nproc) as pool:
            for row in pool.imap_unordered(process_case, tasks):
                record_progress(row)

    logger.info("ALL DONE.")
    logger.info("OK: %d", ok_count)
    logger.info("Skipped: %d", skipped_count)
    logger.info("Failed: %d", failed_count)
    logger.info("Manifest: %s", manifest_path)
    logger.info("Log file: %s", log_file)
    if failed_count and args.fail_on_error:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
