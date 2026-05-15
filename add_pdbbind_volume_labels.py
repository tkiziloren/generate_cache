import argparse
import csv
import os
import shutil
from dataclasses import dataclass
from typing import Iterable

import h5py
import numpy as np


VDW_RADII = {
    "H": 1.20,
    "B": 1.92,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "CL": 1.75,
    "BR": 1.85,
    "I": 1.98,
    "SE": 1.90,
}

DEFAULT_RADIUS = 1.70


@dataclass(frozen=True)
class AtomSet:
    coords: np.ndarray
    elements: list[str]


@dataclass(frozen=True)
class GridSpec:
    origin: np.ndarray
    box_size: int
    resolution: float

    @property
    def shape(self):
        return (self.box_size, self.box_size, self.box_size)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Add ligand-volume/cavity-volume labels to existing PDBBind gridfix H5 caches. "
            "The original cache files are copied by default; use --in-place only for deliberate schema updates."
        )
    )
    parser.add_argument("--dataset-root", default="/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set/")
    parser.add_argument(
        "--h5-dir",
        default="/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set_minimal_cache_gridfix_v1/box72",
    )
    parser.add_argument(
        "--output-dir",
        default="/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set_minimal_cache_gridfix_v1/box72_volume_labels",
    )
    parser.add_argument("--case-list", default=None)
    parser.add_argument("--cases", nargs="*", default=None)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--in-place", action="store_true")
    parser.add_argument("--overwrite-labels", action="store_true")
    parser.add_argument("--include-hydrogen", action="store_true")
    parser.add_argument(
        "--ligand-padding",
        type=float,
        default=1.5,
        help="Extra Angstrom radius added to each ligand heavy atom for the ligand-volume label.",
    )
    parser.add_argument(
        "--protein-padding",
        type=float,
        default=0.0,
        help="Extra Angstrom radius added to each protein heavy atom when excluding protein volume.",
    )
    parser.add_argument("--ligand-label-name", default="binding_site_ligand_volume")
    parser.add_argument("--cavity-label-name", default="binding_site_cavity_volume")
    parser.add_argument("--summary-csv", default=None)
    return parser.parse_args()


def list_cases(h5_dir: str, case_list: str | None, explicit_cases: Iterable[str] | None, limit: int | None):
    if explicit_cases:
        cases = [case[:-3] if case.endswith(".h5") else case for case in explicit_cases]
    elif case_list:
        with open(case_list, "r") as handle:
            cases = [line.strip() for line in handle if line.strip()]
    else:
        cases = sorted(os.path.splitext(name)[0] for name in os.listdir(h5_dir) if name.endswith(".h5"))

    if limit is not None:
        cases = cases[:limit]
    return cases


def normalize_element(raw: str) -> str:
    token = "".join(ch for ch in raw.strip() if ch.isalpha())
    if not token:
        return ""
    token = token.upper()
    if len(token) >= 2 and token[:2] in VDW_RADII:
        return token[:2]
    return token[0]


def parse_mol2_atoms(path: str, include_hydrogen: bool = False) -> AtomSet:
    coords = []
    elements = []
    in_atoms = False

    with open(path, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                continue
            if line.startswith("@<TRIPOS>") and in_atoms:
                break
            if not in_atoms:
                continue

            parts = line.split()
            if len(parts) < 6:
                continue
            element = normalize_element(parts[5].split(".")[0] or parts[1])
            if not include_hydrogen and element == "H":
                continue
            coords.append([float(parts[2]), float(parts[3]), float(parts[4])])
            elements.append(element)

    if not coords:
        raise RuntimeError(f"No ligand atoms parsed from {path}")
    return AtomSet(coords=np.asarray(coords, dtype=np.float64), elements=elements)


def parse_pdb_atoms(path: str, include_hydrogen: bool = False) -> AtomSet:
    coords = []
    elements = []

    with open(path, "r") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if len(line) < 54:
                continue

            element = normalize_element(line[76:78] if len(line) >= 78 else "")
            if not element:
                element = normalize_element(line[12:16])
            if not include_hydrogen and element == "H":
                continue

            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            elements.append(element)

    if not coords:
        raise RuntimeError(f"No atoms parsed from {path}")
    return AtomSet(coords=np.asarray(coords, dtype=np.float64), elements=elements)


def read_grid_spec(h5f: h5py.File) -> GridSpec:
    box_size = int(h5f.attrs["box_size"])
    resolution = float(h5f.attrs["resolution"])
    origin = np.asarray(h5f.attrs["grid_origin"], dtype=np.float64)
    return GridSpec(origin=origin, box_size=box_size, resolution=resolution)


def atoms_to_volume_mask(atom_set: AtomSet, grid: GridSpec, padding: float) -> np.ndarray:
    mask = np.zeros(grid.shape, dtype=bool)
    axes = [grid.origin[axis] + np.arange(grid.box_size) * grid.resolution for axis in range(3)]

    for coord, element in zip(atom_set.coords, atom_set.elements):
        radius = VDW_RADII.get(element, DEFAULT_RADIUS) + padding
        lo = np.floor((coord - radius - grid.origin) / grid.resolution).astype(int)
        hi = np.ceil((coord + radius - grid.origin) / grid.resolution).astype(int) + 1
        lo = np.maximum(lo, 0)
        hi = np.minimum(hi, grid.box_size)

        if np.any(lo >= hi):
            continue

        dx2 = (axes[0][lo[0] : hi[0]] - coord[0]) ** 2
        dy2 = (axes[1][lo[1] : hi[1]] - coord[1]) ** 2
        dz2 = (axes[2][lo[2] : hi[2]] - coord[2]) ** 2
        local = dx2[:, None, None] + dy2[None, :, None] + dz2[None, None, :] <= radius * radius
        mask[lo[0] : hi[0], lo[1] : hi[1], lo[2] : hi[2]] |= local

    return mask


def dice(a: np.ndarray, b: np.ndarray) -> float:
    a_sum = int(a.sum())
    b_sum = int(b.sum())
    if a_sum + b_sum == 0:
        return 1.0
    return float(2 * np.logical_and(a, b).sum() / (a_sum + b_sum))


def jaccard(a: np.ndarray, b: np.ndarray) -> float:
    union = int(np.logical_or(a, b).sum())
    if union == 0:
        return 1.0
    return float(np.logical_and(a, b).sum() / union)


def write_or_replace_label(h5f: h5py.File, name: str, data: np.ndarray, overwrite: bool):
    group = h5f.require_group("label")
    if name in group:
        if not overwrite:
            return False
        del group[name]
    group.create_dataset(name, data=data.astype(np.uint8), compression="gzip")
    labels = set()
    if "labels" in h5f.attrs:
        labels.update(str(label) for label in h5f.attrs["labels"])
    labels.add(name)
    h5f.attrs["labels"] = sorted(labels)
    return True


def prepare_output_file(src_h5: str, dst_h5: str, in_place: bool):
    if in_place:
        return src_h5
    os.makedirs(os.path.dirname(dst_h5), exist_ok=True)
    tmp_h5 = dst_h5 + ".tmp"
    if os.path.exists(tmp_h5):
        os.remove(tmp_h5)
    shutil.copy2(src_h5, tmp_h5)
    os.replace(tmp_h5, dst_h5)
    return dst_h5


def process_case(args, case: str):
    src_h5 = os.path.join(args.h5_dir, f"{case}.h5")
    if not os.path.exists(src_h5):
        raise RuntimeError(f"Missing H5 file: {src_h5}")

    dst_h5 = src_h5 if args.in_place else os.path.join(args.output_dir, f"{case}.h5")
    path = prepare_output_file(src_h5, dst_h5, args.in_place)

    ligand_mol2 = os.path.join(args.dataset_root, case, f"{case}_ligand.mol2")
    protein_pdb = os.path.join(args.dataset_root, case, f"{case}_protein.pdb")
    if not os.path.exists(ligand_mol2):
        raise RuntimeError(f"Missing ligand mol2: {ligand_mol2}")
    if not os.path.exists(protein_pdb):
        raise RuntimeError(f"Missing protein pdb: {protein_pdb}")

    ligand_atoms = parse_mol2_atoms(ligand_mol2, include_hydrogen=args.include_hydrogen)
    protein_atoms = parse_pdb_atoms(protein_pdb, include_hydrogen=args.include_hydrogen)

    with h5py.File(path, "a") as h5f:
        grid = read_grid_spec(h5f)
        ligand_volume = atoms_to_volume_mask(ligand_atoms, grid, padding=args.ligand_padding)
        protein_volume = atoms_to_volume_mask(protein_atoms, grid, padding=args.protein_padding)
        cavity_volume = np.logical_and(ligand_volume, ~protein_volume)

        wrote_ligand = write_or_replace_label(
            h5f, args.ligand_label_name, ligand_volume, overwrite=args.overwrite_labels
        )
        wrote_cavity = write_or_replace_label(
            h5f, args.cavity_label_name, cavity_volume, overwrite=args.overwrite_labels
        )

        h5f.attrs["volume_label_version"] = "pdbbind_ligand_volume_v1"
        h5f.attrs["volume_label_ligand_padding_angstrom"] = args.ligand_padding
        h5f.attrs["volume_label_protein_padding_angstrom"] = args.protein_padding
        h5f.attrs["volume_label_hydrogen_included"] = bool(args.include_hydrogen)
        h5f.attrs["volume_label_ligand_label_name"] = args.ligand_label_name
        h5f.attrs["volume_label_cavity_label_name"] = args.cavity_label_name

        row = {
            "case": case,
            "path": path,
            "box_size": grid.box_size,
            "resolution": grid.resolution,
            "ligand_atoms": len(ligand_atoms.elements),
            "protein_atoms": len(protein_atoms.elements),
            "ligand_volume_voxels": int(ligand_volume.sum()),
            "cavity_volume_voxels": int(cavity_volume.sum()),
            "protein_overlap_removed_voxels": int(np.logical_and(ligand_volume, protein_volume).sum()),
            "wrote_ligand_label": wrote_ligand,
            "wrote_cavity_label": wrote_cavity,
        }

        for existing_name in ("binding_site_calculated", "binding_site_in_dataset"):
            if f"label/{existing_name}" in h5f:
                existing = h5f[f"label/{existing_name}"][:].astype(bool)
                row[f"{existing_name}_voxels"] = int(existing.sum())
                row[f"{existing_name}_dice_cavity"] = dice(existing, cavity_volume)
                row[f"{existing_name}_iou_cavity"] = jaccard(existing, cavity_volume)

        return row


def main():
    args = parse_args()
    if args.in_place:
        args.output_dir = args.h5_dir

    cases = list_cases(args.h5_dir, args.case_list, args.cases, args.limit)
    if not cases:
        raise SystemExit(f"No cases found in {args.h5_dir}")

    rows = []
    failures = []
    for case in cases:
        try:
            row = process_case(args, case)
            rows.append(row)
            print(
                f"[OK] {case}: ligand={row['ligand_volume_voxels']} "
                f"cavity={row['cavity_volume_voxels']} "
                f"removed={row['protein_overlap_removed_voxels']}",
                flush=True,
            )
        except Exception as exc:
            failures.append((case, str(exc)))
            print(f"[FAIL] {case}: {exc}", flush=True)

    if rows:
        csv_path = args.summary_csv or os.path.join(args.output_dir, "volume_label_summary.csv")
        os.makedirs(os.path.dirname(csv_path), exist_ok=True)
        fieldnames = sorted({key for row in rows for key in row})
        with open(csv_path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print(f"Summary CSV written: {csv_path}")

    if failures:
        fail_path = os.path.join(args.output_dir, "volume_label_failures.txt")
        os.makedirs(os.path.dirname(fail_path), exist_ok=True)
        with open(fail_path, "w") as handle:
            for case, error in failures:
                handle.write(f"{case}\t{error}\n")
        raise SystemExit(f"Failed cases: {len(failures)}; see {fail_path}")


if __name__ == "__main__":
    main()
