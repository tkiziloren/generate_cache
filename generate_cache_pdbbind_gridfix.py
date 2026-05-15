import argparse
import contextlib
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from multiprocessing import Pool
from typing import Iterable, Optional

import h5py
import numpy as np
import prody
import prody as pr
from scipy.ndimage import map_coordinates
from scipy.spatial import cKDTree

import openbabel.pybel as pybel
from tfbio.data import Featurizer


BOX_SIZE_DEFAULT = 72
APBS_GRID_POINTS = 161
APBS_SPAN_ANGSTROM = 160.0
GRIDFIX_VERSION = "gridfix_v1"

HYDROPHOBICITY = {
    "ALA": 1.8, "CYS": 2.5, "ASP": -3.5, "GLU": -3.5, "PHE": 2.8,
    "GLY": -0.4, "HIS": -3.2, "ILE": 4.5, "LYS": -3.9, "LEU": 3.8,
    "MET": 1.9, "ASN": -3.5, "PRO": -1.6, "GLN": -3.5, "ARG": -4.5,
    "SER": -0.8, "THR": -0.7, "VAL": 4.2, "TRP": -0.9, "TYR": -1.3,
}

STANDARD_PROTEIN_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "ASH", "GLH", "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CYX", "CYM",
}


@contextlib.contextmanager
def suppress_native_stderr():
    original_stderr = os.dup(2)
    try:
        with open(os.devnull, "w") as devnull:
            os.dup2(devnull.fileno(), 2)
            yield
    finally:
        os.dup2(original_stderr, 2)
        os.close(original_stderr)


@dataclass(frozen=True)
class GridSpec:
    center: np.ndarray
    box_size: int
    resolution: float
    origin: np.ndarray

    @property
    def shape(self):
        return (self.box_size, self.box_size, self.box_size)

    @property
    def max_point(self):
        return self.origin + (self.box_size - 1) * self.resolution

    @property
    def span(self):
        return (self.box_size - 1) * self.resolution


@dataclass(frozen=True)
class DxGrid:
    data: np.ndarray
    origin: np.ndarray
    deltas: np.ndarray

    @property
    def resolution(self):
        nonzero = np.abs(self.deltas).sum(axis=1)
        if not np.allclose(nonzero, nonzero[0]):
            raise ValueError(f"Non-uniform DX deltas are not supported: {self.deltas}")
        return float(nonzero[0])


def make_grid_spec(center: np.ndarray, box_size: int, resolution: Optional[float] = None) -> GridSpec:
    if box_size < 2:
        raise ValueError("box_size must be at least 2")
    if resolution is None:
        resolution = APBS_SPAN_ANGSTROM / float(box_size - 1)
    center = np.asarray(center, dtype=np.float64)
    origin = center - 0.5 * (box_size - 1) * float(resolution)
    return GridSpec(center=center, box_size=box_size, resolution=float(resolution), origin=origin)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate PDBBind HDF5 cache with fixed point-grid semantics. "
            "APBS is computed on 161 grid points spanning 160 Angstrom and "
            "resampled onto the requested target grid."
        )
    )
    parser.add_argument("--dataset-root", default="/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set/")
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--case-list", default="pdbbind_gridfix160_full_fits.txt")
    parser.add_argument("--cases", nargs="*", default=None, help="Optional explicit case ids, e.g. 1a4w 1a1e")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--box-size", type=int, default=BOX_SIZE_DEFAULT)
    parser.add_argument("--resolution", type=float, default=None)
    parser.add_argument("--nproc", type=int, default=1)
    parser.add_argument("--max-size", type=int, default=None, help="Only used to derive case-list if --case-list is omitted")
    parser.add_argument("--apbs-bin", default="apbs")
    parser.add_argument("--obabel-bin", default="obabel")
    parser.add_argument("--pdb2pqr-bin", default="pdb2pqr30")
    parser.add_argument("--apbs-timeout", type=int, default=120)
    parser.add_argument("--keep-temp", action="store_true")
    parser.add_argument("--fail-on-error", action="store_true")
    return parser.parse_args()


def default_output_dir(box_size: int) -> str:
    return (
        "/Users/tevfik/Sandbox/github/PHD/data/pdbbind/"
        f"refined-set_minimal_cache_gridfix_v1/box{box_size}/"
    )


def read_cases(case_list: str, explicit_cases: Optional[Iterable[str]], limit: Optional[int]):
    if explicit_cases:
        cases = [case.strip() for case in explicit_cases if case.strip()]
    else:
        with open(case_list, "r") as handle:
            cases = [line.strip() for line in handle if line.strip()]
    if limit is not None:
        cases = cases[:limit]
    return cases


def mol2_to_pdb(mol2_file, pdb_file, obabel_bin="obabel", timeout=60):
    cmd = [obabel_bin, "-imol2", mol2_file, "-opdb", "-O", pdb_file]
    try:
        result = subprocess.run(cmd, capture_output=True, timeout=timeout, text=True)
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f"obabel conversion timed out for {mol2_file}") from exc
    if result.returncode != 0:
        raise RuntimeError(f"obabel conversion failed for {mol2_file}: {result.stderr}")
    if not os.path.exists(pdb_file) or os.path.getsize(pdb_file) == 0:
        raise RuntimeError(f"obabel created empty or no file: {pdb_file}")


def select_protein_heavy_atoms(atom_group, source_name):
    resnames = np.char.upper(np.asarray(atom_group.getResnames(), dtype=str))
    elements = np.char.upper(np.asarray(atom_group.getElements(), dtype=str))
    names = np.asarray(atom_group.getNames(), dtype=str)

    is_standard_residue = np.isin(resnames, list(STANDARD_PROTEIN_RESNAMES))
    is_hydrogen_name = np.asarray(
        [
            name.upper().startswith("H")
            or (len(name) > 1 and name[0].isdigit() and name[1].upper() == "H")
            for name in names
        ],
        dtype=bool,
    )
    is_hydrogen = (elements == "H") | is_hydrogen_name
    indices = np.where(is_standard_residue & ~is_hydrogen)[0]
    if indices.size == 0:
        raise RuntimeError(f"No standard protein heavy atoms found in: {source_name}")
    return atom_group[indices]


def pdb2pqr_wrapper(selected_pdb, pqr_file, pdb2pqr_bin="pdb2pqr30"):
    cmd = [
        pdb2pqr_bin,
        "--with-ph=7.4",
        "--ff=AMBER",
        "--keep-chain",
        selected_pdb,
        pqr_file,
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"pdb2pqr failed for {selected_pdb}: {err.decode(errors='replace')}")
    if not os.path.exists(pqr_file) or os.path.getsize(pqr_file) == 0:
        raise RuntimeError(f"pdb2pqr created empty or no file: {pqr_file}")
    return pqr_file


def prepare_protein_pocket(protein_pdb, ligand_pdb, pocket_pdb, selected_pdb):
    structure = pr.parsePDB(protein_pdb)
    if structure is None or structure.numAtoms() == 0:
        raise RuntimeError(f"No atoms found in protein PDB: {protein_pdb}")

    structure = select_protein_heavy_atoms(structure, protein_pdb)
    pr.writePDB(protein_pdb, structure)

    protein = pr.parsePDB(protein_pdb)
    ligand = pr.parsePDB(ligand_pdb)
    if ligand is None or ligand.numAtoms() == 0:
        raise RuntimeError(f"No ligand atoms found in: {ligand_pdb}")

    lresname = ligand.getResnames()[0]
    complex_atoms = ligand + protein
    selected_chains = complex_atoms.select(f"same chain as exwithin 7 of resname {lresname}")
    if selected_chains is None or selected_chains.numAtoms() == 0:
        raise RuntimeError(f"No selected chains near ligand for {protein_pdb}")

    selected = select_protein_heavy_atoms(selected_chains, protein_pdb)
    pr.writePDB(selected_pdb, selected)

    with open(selected_pdb, "r") as fin:
        selected_lines = fin.read().splitlines(True)
    with open(selected_pdb, "w") as fout:
        fout.writelines(selected_lines[1:])

    selected = pr.parsePDB(selected_pdb)
    complex_atoms = ligand + selected
    pocket = complex_atoms.select(f"same residue as exwithin 4.5 of resname {lresname}")
    if pocket is None or pocket.numAtoms() == 0:
        raise RuntimeError(f"No calculated pocket atoms for {protein_pdb}")
    pr.writePDB(pocket_pdb, pocket)
    return selected_pdb, pocket_pdb


def protein_center(protein_pdb):
    atoms = prody.parsePDB(protein_pdb).getCoords()
    return atoms.mean(axis=0)


def load_dx_grid(dx_file) -> DxGrid:
    counts = None
    origin = None
    deltas = []
    values = []
    items = None
    reading_data = False

    with open(dx_file, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 8 and parts[:4] == ["object", "1", "class", "gridpositions"]:
                counts = tuple(int(x) for x in parts[-3:])
                continue
            if parts[0] == "origin":
                origin = np.array([float(x) for x in parts[1:4]], dtype=np.float64)
                continue
            if parts[0] == "delta":
                deltas.append([float(x) for x in parts[1:4]])
                continue
            if "data follows" in line:
                reading_data = True
                if "items" in parts:
                    items = int(parts[parts.index("items") + 1])
                continue
            if reading_data:
                try:
                    values.extend(float(x) for x in parts)
                except ValueError:
                    continue
                if items is not None and len(values) >= items:
                    break

    if counts is None or origin is None or len(deltas) != 3:
        raise ValueError(f"Could not parse DX grid metadata from {dx_file}")

    expected = int(np.prod(counts))
    if len(values) < expected:
        raise ValueError(f"DX grid has {len(values)} values, expected {expected}: {dx_file}")

    data = np.array(values[:expected], dtype=np.float32).reshape(counts)
    return DxGrid(data=data, origin=origin, deltas=np.asarray(deltas, dtype=np.float64))


def resample_dx_to_grid(dx_grid: DxGrid, target: GridSpec) -> np.ndarray:
    source_resolution = dx_grid.resolution
    if dx_grid.data.shape == target.shape and np.allclose(dx_grid.origin, target.origin):
        return dx_grid.data.astype(np.float32, copy=True)

    axes = []
    for axis in range(3):
        target_coords = target.origin[axis] + np.arange(target.box_size) * target.resolution
        source_index = (target_coords - dx_grid.origin[axis]) / source_resolution
        axes.append(source_index)
    mesh = np.meshgrid(*axes, indexing="ij")
    resampled = map_coordinates(dx_grid.data, mesh, order=1, mode="nearest")
    return resampled.astype(np.float32)


def run_apbs(protein_pdb, protein_pqr, work_dir, case_name, target: GridSpec, apbs_bin="apbs", timeout=120):
    apbs_case_dir = os.path.join(work_dir, "apbs")
    os.makedirs(apbs_case_dir, exist_ok=True)
    shutil.copy2(protein_pqr, apbs_case_dir)

    apbs_in_file = os.path.join(apbs_case_dir, "apbs.in")
    apbs_out_file = os.path.join(apbs_case_dir, "apbs_grid.dx")
    center = protein_center(protein_pdb)

    with open(apbs_in_file, "w") as handle:
        handle.write(f"""read
    mol pqr {os.path.basename(protein_pqr)}
end
elec
    mg-manual
    mol 1
    dime {APBS_GRID_POINTS} {APBS_GRID_POINTS} {APBS_GRID_POINTS}
    glen {APBS_SPAN_ANGSTROM:.6f} {APBS_SPAN_ANGSTROM:.6f} {APBS_SPAN_ANGSTROM:.6f}
    gcent {center[0]:.6f} {center[1]:.6f} {center[2]:.6f}
    mol 1
    lpbe
    bcfl mdh
    ion charge 1 conc 0.100 radius 2.0
    ion charge -1 conc 0.100 radius 2.0
    pdie 4.0
    sdie 78.54
    srfm smol
    srad 1.4
    sdens 10.0
    chgm spl2
    swin 0.3
    temp 298.15
    calcenergy total
    calcforce no
    write pot dx apbs_grid
end
quit
""")

    try:
        result = subprocess.run(
            [apbs_bin, apbs_in_file],
            cwd=apbs_case_dir,
            check=True,
            timeout=timeout,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"APBS failed for {case_name} with exit code {exc.returncode}\n"
            f"input={apbs_in_file}\nstdout={exc.stdout}\nstderr={exc.stderr}"
        ) from exc
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f"APBS timed out for {case_name} after {timeout}s") from exc

    dx_grid = load_dx_grid(apbs_out_file)
    return resample_dx_to_grid(dx_grid, target)


def coords_to_indices(coords: np.ndarray, grid: GridSpec):
    indices = np.rint((coords - grid.origin) / grid.resolution).astype(np.int64)
    in_box = ((indices >= 0) & (indices < grid.box_size)).all(axis=1)
    return indices, in_box


def pdb_to_mask(pdb_file, grid: GridSpec) -> np.ndarray:
    mol = prody.parsePDB(pdb_file)
    if mol is None or mol.numAtoms() == 0:
        raise RuntimeError(f"No atoms found in {pdb_file}")
    coords = mol.getCoords()
    indices, in_box = coords_to_indices(coords, grid)
    mask = np.zeros(grid.shape, dtype=np.uint8)
    for x, y, z in indices[in_box]:
        mask[x, y, z] = 1
    return mask


def compute_distance_grid(atom_pdb, grid: GridSpec) -> np.ndarray:
    mol = prody.parsePDB(atom_pdb)
    if mol is None or mol.numAtoms() == 0:
        raise RuntimeError(f"No atoms found in {atom_pdb}")
    coords = mol.getCoords()
    axes = [grid.origin[axis] + np.arange(grid.box_size) * grid.resolution for axis in range(3)]
    mesh = np.meshgrid(*axes, indexing="ij")
    grid_points = np.stack([axis.ravel() for axis in mesh], axis=1)
    dists, _ = cKDTree(coords).query(grid_points)
    return dists.reshape(grid.shape).astype(np.float32)


def compute_hydrophobicity_grid(protein_pdb, grid: GridSpec) -> np.ndarray:
    mol = prody.parsePDB(protein_pdb)
    if mol is None or mol.numAtoms() == 0:
        raise RuntimeError(f"No atoms found in {protein_pdb}")
    coords = mol.getCoords()
    residues = mol.getResnames()
    indices, in_box = coords_to_indices(coords, grid)
    hydro_grid = np.zeros(grid.shape, dtype=np.float32)
    for (x, y, z), res in zip(indices[in_box], residues[in_box]):
        hydro_grid[x, y, z] = HYDROPHOBICITY.get(res, 0.0)
    return hydro_grid


def featurize_protein_on_grid(protein_pdb, grid: GridSpec, featurizer: Featurizer) -> np.ndarray:
    with suppress_native_stderr():
        molecule = next(pybel.readfile("pdb", protein_pdb))
    coords, features = featurizer.get_features(molecule, molcode=-1.0)
    indices, in_box = coords_to_indices(coords, grid)

    out = np.zeros((features.shape[1],) + grid.shape, dtype=np.float32)
    for (x, y, z), feature_vector in zip(indices[in_box], features[in_box]):
        out[:, x, y, z] += feature_vector
    return out


def save_atomic_features_separate(h5f, atomic_features, feature_names):
    features_group = h5f.require_group("features")
    expected = atomic_features.shape[1:]
    for idx, name in enumerate(feature_names):
        channel = atomic_features[idx]
        if channel.shape != expected:
            raise RuntimeError(f"Atomic channel shape mismatch for {name}: {channel.shape} vs {expected}")
        features_group.create_dataset(f"atomic_{name}", data=channel, compression="gzip")


def assert_cache_shapes(feature_arrays, label_arrays, grid: GridSpec):
    for name, array in feature_arrays.items():
        if array.shape != grid.shape:
            raise RuntimeError(f"Feature {name} shape {array.shape} does not match target {grid.shape}")
    for name, array in label_arrays.items():
        if array.shape != grid.shape:
            raise RuntimeError(f"Label {name} shape {array.shape} does not match target {grid.shape}")


def process_pdbbind_case(
    case_dir,
    output_h5,
    output_dir,
    case_name,
    box_size,
    resolution,
    apbs_bin,
    obabel_bin,
    pdb2pqr_bin,
    apbs_timeout,
    keep_temp,
):
    if os.path.exists(output_h5):
        print(f"[SKIP] Already exists: {output_h5}", flush=True)
        return

    tmp_root = os.path.join(output_dir, "_tmp")
    os.makedirs(tmp_root, exist_ok=True)
    work_dir = tempfile.mkdtemp(prefix=f"{case_name}_", dir=tmp_root)
    tmp_h5 = output_h5 + ".tmp"

    try:
        src_protein_pdb = os.path.join(case_dir, f"{case_name}_protein.pdb")
        src_ligand_mol2 = os.path.join(case_dir, f"{case_name}_ligand.mol2")
        src_site_pdb = os.path.join(case_dir, f"{case_name}_pocket.pdb")

        protein_pdb = os.path.join(work_dir, f"{case_name}_protein.pdb")
        ligand_pdb = os.path.join(work_dir, f"{case_name}_ligand.pdb")
        selected_pdb = os.path.join(work_dir, f"{case_name}_selected.pdb")
        protein_pqr = os.path.join(work_dir, f"{case_name}_selected.pqr")
        calculated_pocket_pdb = os.path.join(work_dir, f"{case_name}_calculated_pocket.pdb")

        shutil.copy2(src_protein_pdb, protein_pdb)
        mol2_to_pdb(src_ligand_mol2, ligand_pdb, obabel_bin=obabel_bin)
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
            "binding_site_in_dataset": pdb_to_mask(src_site_pdb, grid),
        }
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

            h5f.attrs["schema_version"] = GRIDFIX_VERSION
            h5f.attrs["grid_convention"] = "point_grid_centered_configurable_span"
            h5f.attrs["box_size"] = box_size
            h5f.attrs["resolution"] = grid.resolution
            h5f.attrs["center"] = grid.center.tolist()
            h5f.attrs["grid_origin"] = grid.origin.tolist()
            h5f.attrs["grid_max"] = grid.max_point.tolist()
            h5f.attrs["physical_span_angstrom"] = grid.span
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

        os.replace(tmp_h5, output_h5)
        print(f"[OK] {case_name} -> {output_h5}", flush=True)
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
        apbs_bin,
        obabel_bin,
        pdb2pqr_bin,
        apbs_timeout,
        keep_temp,
        failed_cases_path,
    ) = task
    try:
        process_pdbbind_case(
            case_dir=case_dir,
            output_h5=output_h5,
            output_dir=output_dir,
            case_name=case_name,
            box_size=box_size,
            resolution=resolution,
            apbs_bin=apbs_bin,
            obabel_bin=obabel_bin,
            pdb2pqr_bin=pdb2pqr_bin,
            apbs_timeout=apbs_timeout,
            keep_temp=keep_temp,
        )
        return True
    except Exception as exc:
        print(f"[FAILED] {case_name}: {exc}", flush=True)
        with open(failed_cases_path, "a") as handle:
            handle.write(f"{case_name}\t{exc}\n")
        if os.path.exists(output_h5):
            os.remove(output_h5)
        return False


def main():
    args = parse_args()
    if args.case_list is None and args.max_size is not None:
        args.case_list = f"pdbbind_{args.max_size}_fits.txt"

    output_dir = args.output_dir or default_output_dir(args.box_size)
    os.makedirs(output_dir, exist_ok=True)

    resolution = args.resolution
    if resolution is None:
        resolution = APBS_SPAN_ANGSTROM / float(args.box_size - 1)

    failed_cases_path = os.path.join(output_dir, "failed_cases.txt")
    if os.path.exists(failed_cases_path):
        os.remove(failed_cases_path)

    cases = read_cases(args.case_list, args.cases, args.limit)
    print(f"Gridfix version: {GRIDFIX_VERSION}")
    print(f"Dataset root: {args.dataset_root}")
    print(f"Output dir: {output_dir}")
    print(f"Cases: {len(cases)}")
    print(f"Box size: {args.box_size}")
    print(f"Resolution: {resolution:.8f} A/voxel")
    print(f"Target point span: {(args.box_size - 1) * resolution:.3f} A")
    print(f"APBS source: {APBS_GRID_POINTS} points over {APBS_SPAN_ANGSTROM:.1f} A")

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
                args.apbs_bin,
                args.obabel_bin,
                args.pdb2pqr_bin,
                args.apbs_timeout,
                args.keep_temp,
                failed_cases_path,
            )
        )

    if args.nproc == 1:
        results = []
        for task in tasks:
            results.append(process_case(task))
    else:
        with Pool(processes=args.nproc) as pool:
            results = list(pool.imap_unordered(process_case, tasks))

    failed_count = sum(1 for ok in results if not ok)
    print("ALL DONE.")
    print(f"Failed count: {failed_count}")
    print(f"Failed cases written to: {failed_cases_path}")
    if failed_count and args.fail_on_error:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
