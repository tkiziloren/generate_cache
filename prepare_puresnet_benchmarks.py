import argparse
import csv
import math
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np


PU_RESNET_ROOT_DEFAULT = Path("/Users/tevfik/Sandbox/github/PHD/data/external_benchmarks/puresnet")
P2RANK_DATASETS_ROOT_DEFAULT = Path("/Users/tevfik/Sandbox/github/PHD/data/external_benchmarks/p2rank-datasets")
OUTPUT_ROOT_DEFAULT = Path("/Users/tevfik/Sandbox/github/PHD/data/external_benchmarks/puresnet_prepared")
OBABEL_DEFAULT = Path("/Users/tevfik/Sandbox/Tevfik/Projects/phd_examples/generate_cache/.conda/bin/obabel")

MANIFEST_FIELDS = [
    "dataset",
    "case_id",
    "p2rank_case_id",
    "protein_pdb",
    "ligand_pdb",
    "puresnet_protein_source",
    "puresnet_ligand_source",
    "p2rank_protein_source",
    "fpocket_output_pdb",
    "selected_pocket_index",
    "selected_pocket_atoms_pdb",
    "selected_pocket_vertices_pqr",
    "selection_method",
    "selected_center_dca_angstrom",
    "selected_min_atom_distance_angstrom",
    "selected_pocket_atom_count",
    "selected_pocket_vertex_count",
    "fpocket_pocket_count",
    "status",
    "error",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Prepare PUResNet Coach420/BU48 benchmark subsets with matching "
            "P2Rank fpocket outputs and selected fpocket pockets."
        )
    )
    parser.add_argument("--puresnet-root", type=Path, default=PU_RESNET_ROOT_DEFAULT)
    parser.add_argument("--p2rank-root", type=Path, default=P2RANK_DATASETS_ROOT_DEFAULT)
    parser.add_argument("--output-root", type=Path, default=OUTPUT_ROOT_DEFAULT)
    parser.add_argument("--obabel-bin", type=Path, default=OBABEL_DEFAULT)
    parser.add_argument(
        "--selection-method",
        choices=("center_dca", "atom_min_distance"),
        default="center_dca",
        help="How to choose the fpocket pocket associated with the ligand.",
    )
    parser.add_argument(
        "--datasets",
        nargs="+",
        choices=("coach420_puresnet", "bu48_puresnet"),
        default=("coach420_puresnet", "bu48_puresnet"),
    )
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def require_file(path: Path, description: str):
    if not path.exists() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty {description}: {path}")


def copy_file(src: Path, dst: Path, overwrite: bool):
    require_file(src, "source file")
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() and not overwrite:
        return
    shutil.copy2(src, dst)


def convert_mol2_to_pdb(src: Path, dst: Path, obabel_bin: Path, overwrite: bool):
    require_file(src, "MOL2 source")
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() and dst.stat().st_size > 0 and not overwrite:
        return
    if not obabel_bin.exists():
        raise RuntimeError(f"Open Babel executable not found: {obabel_bin}")
    cmd = [str(obabel_bin), "-imol2", str(src), "-opdb", "-O", str(dst)]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    require_file(dst, "converted PDB")


def read_pdb_coords(path: Path) -> np.ndarray:
    coords = []
    with path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                coords.append(parse_coordinate_tail(line, path))
    if not coords:
        raise RuntimeError(f"No atom coordinates found in PDB: {path}")
    return np.asarray(coords, dtype=np.float64)


def read_pqr_coords(path: Path) -> np.ndarray:
    coords = []
    with path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            coords.append(parse_coordinate_tail(line, path))
    if not coords:
        raise RuntimeError(f"No atom coordinates found in PQR: {path}")
    return np.asarray(coords, dtype=np.float64)


def parse_coordinate_tail(line: str, path: Path):
    values = re.findall(r"[+-]?\d+\.\d+", line[30:])
    if len(values) < 3:
        raise RuntimeError(f"Cannot parse coordinates in {path}: {line.rstrip()}")
    return float(values[0]), float(values[1]), float(values[2])


def center(coords: np.ndarray) -> np.ndarray:
    return coords.mean(axis=0)


def min_pair_distance(a: np.ndarray, b: np.ndarray) -> float:
    best = math.inf
    chunk_size = 2048
    for start in range(0, len(a), chunk_size):
        chunk = a[start : start + chunk_size]
        dist = np.linalg.norm(chunk[:, None, :] - b[None, :, :], axis=2)
        best = min(best, float(dist.min()))
    return best


def pocket_index(path: Path) -> int:
    match = re.search(r"pocket(\d+)_atm\.pdb$", path.name)
    if not match:
        raise RuntimeError(f"Cannot parse fpocket pocket index from: {path}")
    return int(match.group(1))


def choose_fpocket_pocket(fpocket_dir: Path, ligand_pdb: Path, selection_method: str):
    pockets_dir = fpocket_dir / "pockets"
    require_file(fpocket_dir / f"{fpocket_dir.name.removesuffix('_out')}_out.pdb", "fpocket output PDB")
    pocket_atom_files = sorted(pockets_dir.glob("pocket*_atm.pdb"), key=pocket_index)
    if not pocket_atom_files:
        raise RuntimeError(f"No fpocket pocket atom files found: {pockets_dir}")

    ligand_coords = read_pdb_coords(ligand_pdb)
    ligand_center = center(ligand_coords)

    scored = []
    for atoms_path in pocket_atom_files:
        idx = pocket_index(atoms_path)
        vertices_path = pockets_dir / f"pocket{idx}_vert.pqr"
        atom_coords = read_pdb_coords(atoms_path)
        center_dca = float(np.linalg.norm(center(atom_coords) - ligand_center))
        atom_min = min_pair_distance(atom_coords, ligand_coords)
        score = center_dca if selection_method == "center_dca" else atom_min
        vertex_count = 0
        if vertices_path.exists() and vertices_path.stat().st_size > 0:
            try:
                vertex_count = int(read_pqr_coords(vertices_path).shape[0])
            except RuntimeError:
                vertex_count = 0
        scored.append(
            {
                "index": idx,
                "atoms_path": atoms_path,
                "vertices_path": vertices_path,
                "center_dca": center_dca,
                "atom_min": atom_min,
                "atom_count": int(atom_coords.shape[0]),
                "vertex_count": vertex_count,
                "score": score,
            }
        )

    scored.sort(key=lambda item: (item["score"], item["index"]))
    return scored[0], len(pocket_atom_files)


def p2rank_coach_case_map(p2rank_root: Path):
    mapping = {}
    for pdb_path in (p2rank_root / "coach420").glob("*.pdb"):
        key = pdb_path.stem[:4].lower()
        if key in mapping:
            raise RuntimeError(f"Ambiguous Coach420 PDB prefix {key}: {mapping[key]} and {pdb_path}")
        mapping[key] = pdb_path.stem
    return mapping


def write_case_list(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for row in rows:
            if row["status"] == "ok":
                handle.write(f"{row['case_id']}\n")


def prepare_coach(args):
    source_root = args.puresnet_root / "coach" / "coach420"
    p2_case_map = p2rank_coach_case_map(args.p2rank_root)
    cases = sorted(path.name for path in source_root.iterdir() if path.is_dir())
    if args.limit:
        cases = cases[: args.limit]

    rows = []
    for case_id in cases:
        row = empty_row("coach420_puresnet", case_id, args.selection_method)
        try:
            p2rank_case_id = p2_case_map[case_id.lower()]
            source_case = source_root / case_id
            p2rank_protein = args.p2rank_root / "coach420" / f"{p2rank_case_id}.pdb"
            fpocket_dir = args.p2rank_root / "coach420" / "fpocket" / f"{p2rank_case_id}_out"
            fpocket_output = fpocket_dir / f"{p2rank_case_id}_out.pdb"
            output_case = args.output_root / "coach420_puresnet" / case_id
            protein_dst = output_case / "protein.pdb"
            ligand_dst = output_case / "ligand.pdb"

            copy_file(p2rank_protein, protein_dst, args.overwrite)
            copy_file(source_case / "ligand.pdb", ligand_dst, args.overwrite)
            selected, pocket_count = choose_fpocket_pocket(fpocket_dir, ligand_dst, args.selection_method)
            selected_atoms_dst = output_case / "fpocket_selected_pocket_atoms.pdb"
            selected_vertices_dst = output_case / "fpocket_selected_pocket_vertices.pqr"
            copy_file(selected["atoms_path"], selected_atoms_dst, args.overwrite)
            if selected["vertices_path"].exists():
                copy_file(selected["vertices_path"], selected_vertices_dst, args.overwrite)

            row.update(
                {
                    "p2rank_case_id": p2rank_case_id,
                    "protein_pdb": str(protein_dst),
                    "ligand_pdb": str(ligand_dst),
                    "puresnet_protein_source": str(source_case / "protein.pdb"),
                    "puresnet_ligand_source": str(source_case / "ligand.pdb"),
                    "p2rank_protein_source": str(p2rank_protein),
                    "fpocket_output_pdb": str(fpocket_output),
                    "selected_pocket_index": selected["index"],
                    "selected_pocket_atoms_pdb": str(selected_atoms_dst),
                    "selected_pocket_vertices_pqr": str(selected_vertices_dst) if selected_vertices_dst.exists() else "",
                    "selected_center_dca_angstrom": f"{selected['center_dca']:.4f}",
                    "selected_min_atom_distance_angstrom": f"{selected['atom_min']:.4f}",
                    "selected_pocket_atom_count": selected["atom_count"],
                    "selected_pocket_vertex_count": selected["vertex_count"],
                    "fpocket_pocket_count": pocket_count,
                    "status": "ok",
                }
            )
        except Exception as exc:
            row["status"] = "failed"
            row["error"] = str(exc)
        rows.append(row)
    return rows


def prepare_bu48(args):
    source_root = args.puresnet_root / "BU48" / "bench2"
    cases = sorted(path.name for path in source_root.iterdir() if path.is_dir() and not path.name.startswith("."))
    if args.limit:
        cases = cases[: args.limit]

    rows = []
    for case_id in cases:
        row = empty_row("bu48_puresnet", case_id, args.selection_method)
        try:
            source_case = source_root / case_id
            p2rank_protein = args.p2rank_root / "joined" / "bu48" / f"{case_id}.pdb"
            fpocket_dir = args.p2rank_root / "joined" / "bu48" / "fpocket" / f"{case_id}_out"
            fpocket_output = fpocket_dir / f"{case_id}_out.pdb"
            output_case = args.output_root / "bu48_puresnet" / case_id
            protein_dst = output_case / "protein.pdb"
            ligand_dst = output_case / "ligand.pdb"
            puresnet_protein_dst = output_case / "protein_puresnet_converted.pdb"

            copy_file(p2rank_protein, protein_dst, args.overwrite)
            convert_mol2_to_pdb(source_case / "ligand.mol2", ligand_dst, args.obabel_bin, args.overwrite)
            convert_mol2_to_pdb(source_case / "protein.mol2", puresnet_protein_dst, args.obabel_bin, args.overwrite)
            selected, pocket_count = choose_fpocket_pocket(fpocket_dir, ligand_dst, args.selection_method)
            selected_atoms_dst = output_case / "fpocket_selected_pocket_atoms.pdb"
            selected_vertices_dst = output_case / "fpocket_selected_pocket_vertices.pqr"
            copy_file(selected["atoms_path"], selected_atoms_dst, args.overwrite)
            if selected["vertices_path"].exists():
                copy_file(selected["vertices_path"], selected_vertices_dst, args.overwrite)

            row.update(
                {
                    "p2rank_case_id": case_id,
                    "protein_pdb": str(protein_dst),
                    "ligand_pdb": str(ligand_dst),
                    "puresnet_protein_source": str(source_case / "protein.mol2"),
                    "puresnet_ligand_source": str(source_case / "ligand.mol2"),
                    "p2rank_protein_source": str(p2rank_protein),
                    "fpocket_output_pdb": str(fpocket_output),
                    "selected_pocket_index": selected["index"],
                    "selected_pocket_atoms_pdb": str(selected_atoms_dst),
                    "selected_pocket_vertices_pqr": str(selected_vertices_dst) if selected_vertices_dst.exists() else "",
                    "selected_center_dca_angstrom": f"{selected['center_dca']:.4f}",
                    "selected_min_atom_distance_angstrom": f"{selected['atom_min']:.4f}",
                    "selected_pocket_atom_count": selected["atom_count"],
                    "selected_pocket_vertex_count": selected["vertex_count"],
                    "fpocket_pocket_count": pocket_count,
                    "status": "ok",
                }
            )
        except Exception as exc:
            row["status"] = "failed"
            row["error"] = str(exc)
        rows.append(row)
    return rows


def empty_row(dataset: str, case_id: str, selection_method: str):
    return {field: "" for field in MANIFEST_FIELDS} | {
        "dataset": dataset,
        "case_id": case_id,
        "selection_method": selection_method,
    }


def write_manifest(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def write_readme(path: Path):
    text = """# PUResNet Benchmark Preparation

This directory stages the PUResNet Coach420 and BU48 benchmark subsets together
with matching fpocket predictions from `rdk/p2rank-datasets`.

Sources:
- PUResNet benchmark zip files provide the protein-ligand pairs used by the
  PUResNet paper.
- `rdk/p2rank-datasets` provides fpocket outputs generated with fpocket v1.0
  default parameters, according to its README.

Prepared subsets:
- `coach420_puresnet`: 298 Coach420 cases present in PUResNet's public zip.
- `bu48_puresnet`: 62 BU48 cases present in PUResNet's public zip.

Per-case files:
- `protein.pdb`: protein structure aligned with the P2Rank/fpocket source.
- `ligand.pdb`: ligand structure from the PUResNet benchmark source.
- `fpocket_selected_pocket_atoms.pdb`: selected fpocket pocket atom file.
- `fpocket_selected_pocket_vertices.pqr`: selected fpocket pocket vertex file.

Pocket selection:
- The default selection method is `center_dca`.
- For every fpocket pocket, the script computes the pocket atom centroid and
  selects the pocket with the smallest distance to the ligand centroid.
- The manifest also records the minimum atom-to-atom distance as a diagnostic.

These files are prepared benchmark inputs. They are not HDF5 training caches yet.
The selected fpocket pocket can be converted into a label/mask in a later cache
generation step.
"""
    path.write_text(text)


def main():
    args = parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)

    all_rows = []
    if "coach420_puresnet" in args.datasets:
        all_rows.extend(prepare_coach(args))
    if "bu48_puresnet" in args.datasets:
        all_rows.extend(prepare_bu48(args))

    write_manifest(args.output_root / "manifest.csv", all_rows)
    for dataset in ("coach420_puresnet", "bu48_puresnet"):
        rows = [row for row in all_rows if row["dataset"] == dataset]
        if rows:
            write_case_list(args.output_root / f"{dataset}_cases.txt", rows)
    write_readme(args.output_root / "README.md")

    ok = sum(row["status"] == "ok" for row in all_rows)
    failed = len(all_rows) - ok
    print(f"Prepared {ok} cases under {args.output_root}")
    if failed:
        print(f"Failed {failed} cases. See manifest.csv for details.")


if __name__ == "__main__":
    main()
