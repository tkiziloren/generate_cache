#!/usr/bin/env python3
"""Clean and append APBS representation channels to existing H5 cache files.

This script does not regenerate APBS, labels, or atomic features. It only adds
channels that can be derived safely from already stored legacy H5 datasets.

By default it also moves leakage/debug channels from ``features/`` to
``auxiliary/`` so they are not offered as model input features.

Full-protein APBS raw values must be produced by a separate source-based APBS
append job. When ``features/electrostatic_grid_v2_full_protein_raw`` is present,
this script adds the derived v2 normalization channels automatically.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import shutil
from pathlib import Path

import h5py
import numpy as np

from feature_schema_metadata import apply_feature_schema_metadata


DEFAULT_SAMPLE_DIR = Path(
    "/Users/tevfik/Sandbox/Tevfik/Projects/phd_examples/codon_h5_samples_2026-05-16"
)


LEGACY_APBS_NAME = "electrostatic_grid"
V1_APBS_PREFIX = "electrostatic_grid_v1_ligand_proximal_chains_7A"
V2_APBS_PREFIX = "electrostatic_grid_v2_full_protein"
AUXILIARY_FEATURES = ("dist_to_ligand", "ligand")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Append H5-derived v2 feature channels. Use --in-place on Codon, "
            "or --output-dir locally to keep source files untouched."
        )
    )
    parser.add_argument("h5_paths", nargs="*", help="One or more H5 files.")
    parser.add_argument("--input-dir", default=None, help="Directory to search for H5 files.")
    parser.add_argument("--glob", default="*.h5", help="Glob used with --input-dir.")
    parser.add_argument("--output-dir", default=None, help="Copy H5 files here before augmenting.")
    parser.add_argument("--in-place", action="store_true", help="Modify input H5 files directly.")
    parser.add_argument("--overwrite", action="store_true", help="Replace derived channels if they already exist.")
    parser.add_argument(
        "--nproc",
        type=int,
        default=1,
        help="Number of parallel H5 files to process. Use 1 for sequential processing.",
    )
    parser.add_argument(
        "--full-signed-scale",
        type=float,
        default=150.0,
        help="Scale used for raw APBS signed representation: value / scale.",
    )
    parser.add_argument(
        "--keep-leakage-features",
        action="store_true",
        help="Keep dist_to_ligand and ligand under features/. By default they are moved to auxiliary/.",
    )
    parser.add_argument(
        "--move-legacy-electrostatic-to-auxiliary",
        action="store_true",
        help=(
            "Deprecated. Prefer --drop-legacy-electrostatic. "
            "Move features/electrostatic_grid to auxiliary/electrostatic_grid after "
            "writing electrostatic_grid_v1_ligand_proximal_chains_7A_raw."
        ),
    )
    parser.add_argument(
        "--drop-legacy-electrostatic",
        action="store_true",
        help="Delete features/electrostatic_grid after writing the explicit v1 APBS raw channel.",
    )
    parser.add_argument(
        "--skip-derived-features",
        action="store_true",
        help="Do not add H5-derived geometry/APBS helper features.",
    )
    parser.add_argument(
        "--summary-csv",
        default=None,
        help="CSV path for augmentation summary. Defaults to <output-dir>/v2_augment_summary.csv or current directory.",
    )
    return parser.parse_args()


def discover_h5_paths(args: argparse.Namespace) -> list[Path]:
    paths = [Path(path).expanduser().resolve() for path in args.h5_paths]
    if args.input_dir:
        paths.extend(sorted(Path(args.input_dir).expanduser().resolve().rglob(args.glob)))
    if not paths and DEFAULT_SAMPLE_DIR.exists():
        paths = sorted(DEFAULT_SAMPLE_DIR.rglob(args.glob))

    unique_paths = []
    seen = set()
    for path in paths:
        if path in seen:
            continue
        seen.add(path)
        if not path.exists():
            raise FileNotFoundError(path)
        if path.suffix != ".h5":
            raise ValueError(f"Expected .h5 file: {path}")
        unique_paths.append(path)
    if not unique_paths:
        raise SystemExit("No H5 files found. Pass paths or use --input-dir.")
    return unique_paths


def prepare_targets(paths: list[Path], args: argparse.Namespace) -> list[tuple[Path, Path]]:
    if args.in_place and args.output_dir:
        raise SystemExit("Use either --in-place or --output-dir, not both.")
    if not args.in_place and not args.output_dir:
        raise SystemExit("Pass --output-dir for copied output, or --in-place to modify inputs directly.")

    if args.in_place:
        return [(path, path) for path in paths]

    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    targets = []
    for source in paths:
        target = output_dir / source.name
        if not target.exists() or args.overwrite:
            shutil.copy2(source, target)
        targets.append((source, target))
    return targets


def clip_minmax(values: np.ndarray, limit: float) -> np.ndarray:
    clipped = np.clip(values, -limit, limit)
    return ((clipped + limit) / (2.0 * limit)).astype(np.float32)


def clip_signed(values: np.ndarray, limit: float) -> np.ndarray:
    return (np.clip(values, -limit, limit) / limit).astype(np.float32)


def full_signed(values: np.ndarray, scale: float) -> np.ndarray:
    if scale <= 0:
        raise ValueError("--full-signed-scale must be positive")
    return (values / scale).astype(np.float32)


def positive(values: np.ndarray, limit: float) -> np.ndarray:
    return (np.clip(values, 0.0, limit) / limit).astype(np.float32)


def negative(values: np.ndarray, limit: float) -> np.ndarray:
    return (np.clip(-values, 0.0, limit) / limit).astype(np.float32)


def robust_minmax(values: np.ndarray, percentile: float = 99.5) -> np.ndarray:
    values = np.asarray(values, dtype=np.float32)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return np.zeros_like(values, dtype=np.float32)
    high = float(np.percentile(finite, percentile))
    if high <= 0.0:
        return np.zeros_like(values, dtype=np.float32)
    return np.clip(values / high, 0.0, 1.0).astype(np.float32)


def gradient_magnitude(values: np.ndarray, resolution: float) -> np.ndarray:
    spacing = float(resolution) if resolution and resolution > 0 else 1.0
    gradients = np.gradient(values.astype(np.float32), spacing, edge_order=1)
    magnitude = np.sqrt(sum(component * component for component in gradients))
    return robust_minmax(magnitude)


def write_feature(h5f: h5py.File, name: str, data: np.ndarray, overwrite: bool, created: list[str]) -> None:
    path = f"features/{name}"
    if path in h5f:
        if not overwrite:
            return
        del h5f[path]
    h5f.create_dataset(path, data=data, compression="gzip")
    created.append(name)


def read_feature(h5f: h5py.File, name: str) -> np.ndarray | None:
    path = f"features/{name}"
    if path in h5f:
        return h5f[path][:].astype(np.float32)
    return None


def move_feature_to_auxiliary(
    h5f: h5py.File,
    name: str,
    overwrite: bool,
    moved: list[str],
    skipped: list[str],
) -> None:
    source_path = f"features/{name}"
    target_path = f"auxiliary/{name}"
    if source_path not in h5f:
        return

    auxiliary_group = h5f.require_group("auxiliary")
    if name in auxiliary_group:
        if overwrite:
            del auxiliary_group[name]
            auxiliary_group.create_dataset(name, data=h5f[source_path][:], compression="gzip")
            moved.append(name)
        else:
            skipped.append(f"{target_path}_already_exists")
    else:
        auxiliary_group.create_dataset(name, data=h5f[source_path][:], compression="gzip")
        moved.append(name)

    del h5f[source_path]


def drop_feature(h5f: h5py.File, name: str, dropped: list[str]) -> None:
    path = f"features/{name}"
    if path in h5f:
        del h5f[path]
        dropped.append(name)


def add_apbs_representations(
    h5f: h5py.File,
    prefix: str,
    raw_values: np.ndarray,
    overwrite: bool,
    full_signed_scale: float,
    created: list[str],
) -> None:
    write_feature(h5f, f"{prefix}_raw", raw_values.astype(np.float32), overwrite, created)
    for limit in (5.0, 10.0, 20.0):
        write_feature(
            h5f,
            f"{prefix}_clip{int(limit)}_minmax",
            clip_minmax(raw_values, limit),
            overwrite,
            created,
        )
    write_feature(
        h5f,
        f"{prefix}_clip20_signed",
        clip_signed(raw_values, 20.0),
        overwrite,
        created,
    )
    write_feature(
        h5f,
        f"{prefix}_full_signed150",
        full_signed(raw_values, full_signed_scale),
        overwrite,
        created,
    )
    write_feature(
        h5f,
        f"{prefix}_clip150_signed",
        clip_signed(raw_values, full_signed_scale),
        overwrite,
        created,
    )
    write_feature(
        h5f,
        f"{prefix}_positive_clip20",
        positive(raw_values, 20.0),
        overwrite,
        created,
    )
    write_feature(
        h5f,
        f"{prefix}_negative_clip20",
        negative(raw_values, 20.0),
        overwrite,
        created,
    )


def add_cache_derived_features(
    h5f: h5py.File,
    overwrite: bool,
    created: list[str],
    full_signed_scale: float = 150.0,
) -> None:
    dist_to_surface = read_feature(h5f, "dist_to_surface")
    hydrophobicity = read_feature(h5f, "hydrophobicity")
    resolution = float(h5f.attrs.get("resolution", 1.0))

    surface_proximity = None
    if dist_to_surface is not None:
        surface_proximity = np.exp(-np.clip(dist_to_surface, 0.0, None) / 3.0).astype(np.float32)
        write_feature(h5f, "protein_proximity_exp3", surface_proximity, overwrite, created)
        write_feature(
            h5f,
            "protein_near_shell_0_3A",
            ((dist_to_surface >= 0.0) & (dist_to_surface <= 3.0)).astype(np.float32),
            overwrite,
            created,
        )
        write_feature(
            h5f,
            "protein_near_shell_3_6A",
            ((dist_to_surface > 3.0) & (dist_to_surface <= 6.0)).astype(np.float32),
            overwrite,
            created,
        )

    if surface_proximity is not None and hydrophobicity is not None:
        write_feature(
            h5f,
            "hydrophobicity_surface_weighted",
            (hydrophobicity * surface_proximity).astype(np.float32),
            overwrite,
            created,
        )

    for prefix in (V1_APBS_PREFIX, V2_APBS_PREFIX):
        raw_values = read_feature(h5f, f"{prefix}_raw")
        if raw_values is None:
            continue
        write_feature(
            h5f,
            f"{prefix}_gradient_magnitude_robust",
            gradient_magnitude(raw_values, resolution),
            overwrite,
            created,
        )
        if surface_proximity is not None:
            write_feature(
                h5f,
                f"{prefix}_clip20_signed_surface_weighted",
                (clip_signed(raw_values, 20.0) * surface_proximity).astype(np.float32),
                overwrite,
                created,
            )
            write_feature(
                h5f,
                f"{prefix}_full_signed150_surface_weighted",
                (full_signed(raw_values, full_signed_scale) * surface_proximity).astype(np.float32),
                overwrite,
                created,
            )


def update_feature_attr(h5f: h5py.File) -> None:
    if "features" not in h5f:
        return
    names = sorted(h5f["features"].keys())
    h5f.attrs["features"] = ",".join(names)


def augment_one(path: Path, args: argparse.Namespace) -> dict[str, str | int]:
    created: list[str] = []
    moved: list[str] = []
    dropped: list[str] = []
    skipped: list[str] = []

    with h5py.File(path, "a") as h5f:
        if "features" not in h5f:
            raise RuntimeError(f"No features group found in {path}")

        legacy_apbs = read_feature(h5f, LEGACY_APBS_NAME)
        if legacy_apbs is not None:
            add_apbs_representations(
                h5f,
                V1_APBS_PREFIX,
                legacy_apbs,
                args.overwrite,
                args.full_signed_scale,
                created,
            )
            if args.drop_legacy_electrostatic:
                drop_feature(h5f, LEGACY_APBS_NAME, dropped)
            elif args.move_legacy_electrostatic_to_auxiliary:
                move_feature_to_auxiliary(h5f, LEGACY_APBS_NAME, args.overwrite, moved, skipped)
        else:
            skipped.append(f"features/{LEGACY_APBS_NAME}")

        v2_raw = read_feature(h5f, f"{V2_APBS_PREFIX}_raw")
        if v2_raw is not None:
            add_apbs_representations(
                h5f,
                V2_APBS_PREFIX,
                v2_raw,
                args.overwrite,
                args.full_signed_scale,
                created,
            )
        else:
            skipped.append(f"features/{V2_APBS_PREFIX}_raw")

        if not args.keep_leakage_features:
            for feature_name in AUXILIARY_FEATURES:
                move_feature_to_auxiliary(h5f, feature_name, args.overwrite, moved, skipped)

        if not args.skip_derived_features:
            add_cache_derived_features(h5f, args.overwrite, created, args.full_signed_scale)

        h5f.attrs["apbs_representation_feature_version"] = "apbs_representations_from_h5_v1"
        h5f.attrs["h5_derived_feature_version"] = "h5_derived_geometry_apbs_v1"
        h5f.attrs["apbs_representation_feature_note"] = (
            "V1 APBS representations are derived from existing electrostatic_grid. "
            "V2 full-protein APBS representations are derived only when "
            "electrostatic_grid_v2_full_protein_raw already exists."
        )
        if f"features/{V2_APBS_PREFIX}_raw" not in h5f:
            h5f.attrs["v2_full_protein_apbs_status"] = "missing_requires_source_apbs_append"
        update_feature_attr(h5f)
        apply_feature_schema_metadata(h5f)

    return {
        "path": str(path),
        "created_count": len(created),
        "created_features": ";".join(created),
        "moved_to_auxiliary": ";".join(moved),
        "dropped_features": ";".join(dropped),
        "skipped_sources": ";".join(skipped),
    }


def summary_path(args: argparse.Namespace) -> Path:
    if args.summary_csv:
        return Path(args.summary_csv).expanduser().resolve()
    if args.output_dir:
        return Path(args.output_dir).expanduser().resolve() / "v2_augment_summary.csv"
    return Path.cwd() / "v2_augment_summary.csv"


def process_target(task: tuple[int, Path, Path, argparse.Namespace]) -> dict[str, str | int]:
    index, source, target, args = task
    row = augment_one(target, args)
    row["index"] = index
    row["source_path"] = str(source)
    return row


def main() -> None:
    args = parse_args()
    if args.nproc < 1:
        raise SystemExit("--nproc must be at least 1")

    source_paths = discover_h5_paths(args)
    targets = prepare_targets(source_paths, args)

    rows = []
    total_targets = len(targets)
    tasks = [
        (index, source, target, args)
        for index, (source, target) in enumerate(targets, start=1)
    ]

    if args.nproc == 1:
        iterator = map(process_target, tasks)
    else:
        executor = ProcessPoolExecutor(max_workers=args.nproc)
        iterator = executor.map(process_target, tasks)

    for completed, row in enumerate(iterator, start=1):
        rows.append(row)
        print(
            f"[PROGRESS] feature schema {completed}/{total_targets} completed | "
            f"remaining={total_targets - completed} | added={row['created_count']} | "
            f"path={row['path']} | skipped={row['skipped_sources']}",
            flush=True,
        )
    if args.nproc != 1:
        executor.shutdown()

    csv_path = summary_path(args)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_path, "w", newline="") as handle:
        fieldnames = [
            "index",
            "source_path",
            "path",
            "created_count",
            "created_features",
            "moved_to_auxiliary",
            "dropped_features",
            "skipped_sources",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Summary written to: {csv_path}", flush=True)


if __name__ == "__main__":
    main()
