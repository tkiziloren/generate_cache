import os
import numpy as np
import h5py
import prody
from scipy.spatial import cKDTree
from scipy.ndimage import zoom
import subprocess
import shutil, os

BOX_SIZE = 64  # Ana grid boyutun (featurizer ve ML modeli için)
APBS_SIZE = BOX_SIZE + 1  # APBS input için (+1 önemli!)

# Kyte-Doolittle hydrophobicity scale
HYDROPHOBICITY = {
    'ALA': 1.8,  'CYS': 2.5,  'ASP': -3.5, 'GLU': -3.5, 'PHE': 2.8,
    'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,  'LYS': -3.9, 'LEU': 3.8,
    'MET': 1.9,  'ASN': -3.5, 'PRO': -1.6, 'GLN': -3.5, 'ARG': -4.5,
    'SER': -0.8, 'THR': -0.7, 'VAL': 4.2,  'TRP': -0.9, 'TYR': -1.3
}

try:
    from tfbio.data import Featurizer
except ImportError:
    raise ImportError("Please install tfbio: pip install git+https://gitlab.com/cheminfIBB/tfbio.git")

def mol2_to_pdb(mol2_file, pdb_file):
    cmd = ["obabel", "-imol2", mol2_file, "-opdb", "-O", pdb_file]
    result = subprocess.run(cmd, capture_output=True)
    if result.returncode != 0:
        print("OpenBabel error:", result.stderr.decode())
        raise RuntimeError(f"obabel conversion failed for {mol2_file}")

def get_grid_params(protein_pdb, box_size, grid_step=1.0):
    atoms = prody.parsePDB(protein_pdb).getCoords()
    center = atoms.mean(axis=0)
    grid_min = center - (box_size // 2) * grid_step
    grid_max = center + (box_size // 2) * grid_step
    grid_shape = (box_size, box_size, box_size)
    return grid_shape, grid_min, grid_max, center

def run_apbs(protein_pdb, cache_dir, case_name, box_size, grid_step=1.0):
    apbs_case_dir = os.path.join(cache_dir, case_name)
    os.makedirs(apbs_case_dir, exist_ok=True)
    pqr_file = os.path.join(apbs_case_dir, "protein.pqr")
    apbs_in_file = os.path.join(apbs_case_dir, "apbs.in")
    apbs_out_file = os.path.join(apbs_case_dir, "apbs_grid.dx")

    grid_shape, grid_min, grid_max, grid_center = get_grid_params(protein_pdb, box_size, grid_step)
    subprocess.run(["pdb2pqr30", "--ff=AMBER", protein_pdb, pqr_file], check=True)

    cglen = grid_max - grid_min
    dime = [max(33, box_size)] * 3  # APBS min 33, burada doğrudan APBS_SIZE veriyoruz
    with open(apbs_in_file, "w") as f:
        f.write(f"""read
    mol pqr {os.path.basename(pqr_file)}
end
elec name prot
    mg-auto
    dime {dime[0]} {dime[1]} {dime[2]}
    cglen {cglen[0]:.1f} {cglen[1]:.1f} {cglen[2]:.1f}
    fglen {cglen[0]:.1f} {cglen[1]:.1f} {cglen[2]:.1f}
    cgcent {grid_center[0]:.2f} {grid_center[1]:.2f} {grid_center[2]:.2f}
    fgcent {grid_center[0]:.2f} {grid_center[1]:.2f} {grid_center[2]:.2f}
    mol 1
    lpbe
    bcfl sdh
    pdie 2.0
    sdie 78.54
    srfm smol
    chgm spl2
    sdens 10.0
    srad 1.4
    swin 0.3
    temp 298.15
    calcenergy total
    calcforce no
    write pot dx apbs_grid
end
quit
""")
    subprocess.run(["apbs", apbs_in_file], cwd=apbs_case_dir, check=True)
    apbs_grid = load_dx_grid(apbs_out_file)
    apbs_grid_resized = resize_grid(apbs_grid, (BOX_SIZE, BOX_SIZE, BOX_SIZE))
    return apbs_grid_resized, grid_center, grid_shape, grid_min, grid_max

def load_dx_grid(dx_file):
    with open(dx_file) as f:
        lines = f.readlines()
    data = []
    for line in lines:
        if line.startswith("object 3"):
            continue
        if line[0].isdigit() or line[0] == '-' or line.strip().startswith('.'):
            data.extend([float(x) for x in line.strip().split()])
    arr = np.array(data)
    return arr

def resize_grid(arr, target_shape):
    orig_len = arr.shape[0]
    possible_shapes = []
    for n1 in range(32, 200):
        if orig_len % n1 != 0: continue
        for n2 in range(32, 200):
            if (orig_len // n1) % n2 != 0: continue
            n3 = orig_len // (n1 * n2)
            if n1 * n2 * n3 == orig_len and n3 < 200:
                possible_shapes.append((n1, n2, n3))
    if not possible_shapes:
        raise ValueError("Could not infer grid shape from .dx data")
    shape = min(possible_shapes, key=lambda s: sum(abs(s[i] - target_shape[i]) for i in range(3)))
    arr3d = arr.reshape(shape)
    zoom_factors = [target_shape[i] / arr3d.shape[i] for i in range(3)]
    arr3d_resized = zoom(arr3d, zoom_factors, order=1)
    return arr3d_resized.astype(np.float32)

def get_protein_shape_mask(protein_pdb, grid_center, box_size, resolution):
    mol = prody.parsePDB(protein_pdb)
    coords = mol.getCoords()
    grid_shape = (box_size, box_size, box_size)
    mask = np.zeros(grid_shape, dtype=np.uint8)
    grid_start = np.array(grid_center) - (box_size // 2) * resolution
    for atom in coords:
        idx = np.round((atom - grid_start) / resolution).astype(int)
        if np.all((idx >= 0) & (idx < box_size)):
            mask[tuple(idx)] = 1
    return mask

def mol2_to_mask(mol2_file, grid_center, box_size, resolution):
    mol = prody.parsePDB(mol2_file)
    coords = mol.getCoords()
    grid_shape = (box_size, box_size, box_size)
    mask = np.zeros(grid_shape, dtype=np.uint8)
    grid_start = np.array(grid_center) - (box_size // 2) * resolution
    for atom in coords:
        idx = np.round((atom - grid_start) / resolution).astype(int)
        if np.all((idx >= 0) & (idx < box_size)):
            mask[tuple(idx)] = 1
    return mask

def compute_distance_to_ligand(ligand_pdb, grid_center, box_size, resolution):
    ligand = prody.parsePDB(ligand_pdb)
    ligand_coords = ligand.getCoords()
    grid = np.indices((box_size, box_size, box_size)).reshape(3, -1).T
    grid_start = np.array(grid_center) - (box_size // 2) * resolution
    grid_points = grid_start + grid * resolution
    tree = cKDTree(ligand_coords)
    dists, _ = tree.query(grid_points)
    dist_grid = dists.reshape((box_size, box_size, box_size)).astype(np.float32)
    return dist_grid

def featurize_protein(protein_pdb, grid_center, box_size, resolution):
    featurizer = Featurizer()
    channels = featurizer.featurize(protein_pdb, box_size, resolution, grid_center)
    return channels.astype(np.float32)

def save_atomic_features_separate(h5f, atomic_features, feature_names):
    features_grp = h5f.require_group("features")
    for idx, name in enumerate(feature_names):
        dset_name = f"atomic_{name}"
        features_grp.create_dataset(dset_name, data=atomic_features[idx], compression="gzip")

### YENİ FEATURE: Hydrophobicity
def compute_hydrophobicity_grid(protein_pdb, grid_center, box_size, resolution):
    mol = prody.parsePDB(protein_pdb)
    coords = mol.getCoords()
    residues = mol.getResnames()
    grid_shape = (box_size, box_size, box_size)
    grid = np.zeros(grid_shape, dtype=np.float32)
    grid_start = np.array(grid_center) - (box_size // 2) * resolution
    for atom, res in zip(coords, residues):
        idx = np.round((atom - grid_start) / resolution).astype(int)
        if np.all((idx >= 0) & (idx < box_size)):
            hydro = HYDROPHOBICITY.get(res, 0.0)
            grid[tuple(idx)] = hydro
    return grid

### YENİ FEATURE: Donor / Acceptor
def compute_donor_acceptor_grids(protein_pdb, grid_center, box_size, resolution):
    mol = prody.parsePDB(protein_pdb)
    coords = mol.getCoords()
    atomnames = mol.getNames()
    grid_donor = np.zeros((box_size, box_size, box_size), dtype=np.uint8)
    grid_acceptor = np.zeros((box_size, box_size, box_size), dtype=np.uint8)
    grid_start = np.array(grid_center) - (box_size // 2) * resolution
    for atom, atom_name in zip(coords, atomnames):
        idx = np.round((atom - grid_start) / resolution).astype(int)
        if np.all((idx >= 0) & (idx < box_size)):
            # Basit kural: N donor, O acceptor
            if atom_name.startswith("N"):
                grid_donor[tuple(idx)] = 1
            if atom_name.startswith("O"):
                grid_acceptor[tuple(idx)] = 1
    return grid_donor, grid_acceptor

def compute_distance_to_surface(protein_pdb, grid_center, box_size, resolution):
    """
    Her grid noktasının en yakın protein atomuna olan uzaklığını döndürür.
    """
    mol = prody.parsePDB(protein_pdb)
    coords = mol.getCoords()
    grid = np.indices((box_size, box_size, box_size)).reshape(3, -1).T  # (N, 3)
    grid_start = np.array(grid_center) - (box_size // 2) * resolution
    grid_points = grid_start + grid * resolution
    tree = cKDTree(coords)
    dists, _ = tree.query(grid_points)
    dist_grid = dists.reshape((box_size, box_size, box_size)).astype(np.float32)
    return dist_grid

def process_protein_case(case_dir, output_h5, box_size, resolution, output_dir=None, case_name=None):
    protein_mol2 = os.path.join(case_dir, "protein.mol2")
    ligand_mol2 = os.path.join(case_dir, "ligand.mol2")
    site_mol2 = os.path.join(case_dir, "site.mol2")

    protein_pdb = os.path.join(case_dir, "protein.pdb")
    ligand_pdb = os.path.join(case_dir, "ligand.pdb")
    site_pdb = os.path.join(case_dir, "site.pdb")
    mol2_to_pdb(protein_mol2, protein_pdb)
    mol2_to_pdb(ligand_mol2, ligand_pdb)
    mol2_to_pdb(site_mol2, site_pdb)

    print(f"{case_dir}: APBS grid calculation...")
    apbs_grid, grid_center, grid_shape, grid_min, grid_max = run_apbs(
        protein_pdb, output_dir, case_name, box_size=APBS_SIZE, grid_step=resolution
    )

    print(f"{case_dir}: Featurizing protein atomic features...")
    featurizer = Featurizer()
    atomic_features = featurizer.featurize(protein_pdb, box_size, resolution, grid_center)
    print(f"{case_dir}: Hydrophobicity grid...")
    hydro_grid = compute_hydrophobicity_grid(protein_pdb, grid_center, box_size, resolution)
    print(f"{case_dir}: Donor/Acceptor grid...")
    donor_grid, acceptor_grid = compute_donor_acceptor_grids(protein_pdb, grid_center, box_size, resolution)
    print(f"{case_dir}: Protein shape mask...")
    shape_mask = get_protein_shape_mask(protein_pdb, grid_center, box_size, resolution)
    print(f"{case_dir}: Ligand distance map...")
    dist2lig = compute_distance_to_ligand(ligand_pdb, grid_center, box_size, resolution)
    print(f"{case_dir}: Binding site mask...")
    label = mol2_to_mask(site_pdb, grid_center, box_size, resolution)
    print(f"{case_dir}: Distance to surface grid...")
    dist2surf = compute_distance_to_surface(protein_pdb, grid_center, box_size, resolution)

    print(f"{case_dir}: Writing h5...")
    with h5py.File(output_h5, "w") as h5f:
        save_atomic_features_separate(h5f, atomic_features, featurizer.FEATURE_NAMES)
        h5f.create_dataset("features/apbs", data=apbs_grid, compression="gzip")
        h5f.create_dataset("features/shape", data=shape_mask, compression="gzip")
        h5f.create_dataset("features/dist2ligand", data=dist2lig, compression="gzip")
        h5f.create_dataset("features/hydrophobicity", data=hydro_grid, compression="gzip")
        h5f.create_dataset("features/donor", data=donor_grid, compression="gzip")
        h5f.create_dataset("features/acceptor", data=acceptor_grid, compression="gzip")
        h5f.create_dataset("features/dist2surface", data=dist2surf, compression="gzip")
        h5f.create_dataset("label/binding_site", data=label, compression="gzip")
        h5f.attrs["box_size"] = box_size
        h5f.attrs["resolution"] = resolution
        h5f.attrs["center"] = grid_center.tolist()
        h5f.attrs["features"] = ",".join(
            [f"atomic_{x}" for x in featurizer.FEATURE_NAMES] +
            ["apbs", "shape", "dist2ligand", "hydrophobicity", "donor", "acceptor", "dist2surface"]  # , "sasa"
        )
        h5f.attrs["label"] = "binding_site"
        h5f.attrs["apbs_grid_shape"] = [box_size, box_size, box_size]
        h5f.attrs["apbs_grid_min"] = grid_min.tolist()
        h5f.attrs["apbs_grid_max"] = grid_max.tolist()
    print(f"{case_dir}: DONE")
    
    try:
        for f in [protein_pdb, ligand_pdb, site_pdb]:
            if os.path.exists(f):
                os.remove(f)
        apbs_case_dir = os.path.join(output_dir, case_name)
        if os.path.exists(apbs_case_dir):
            shutil.rmtree(apbs_case_dir)  
        print(f"{case_dir}: Temporary files are deleted.")
    except Exception as cleanup_err:
        print(f"Warning: For {case_dir} temporary files couldn't be deleted. Error message: {cleanup_err}")

if __name__ == "__main__":
    import tqdm
    scpdb_dir = "/Users/tevfik/Sandbox/github/PHD/data/scPDB"
    output_dir = "/Users/tevfik/Sandbox/github/PHD/data/scPDB_minimal_cache"
    os.makedirs(output_dir, exist_ok=True)
    all_cases = sorted([f for f in os.listdir(scpdb_dir) if os.path.isdir(os.path.join(scpdb_dir, f))])
    for case in tqdm.tqdm(all_cases):
        case_dir = os.path.join(scpdb_dir, case)
        out_h5 = os.path.join(output_dir, f"{case}.h5")
        try:
            process_protein_case(
                case_dir, out_h5, box_size=BOX_SIZE, resolution=1.0, output_dir=output_dir, case_name=case
            )
        except Exception as e:
            print(f"Top-level error in {case}: {e}")