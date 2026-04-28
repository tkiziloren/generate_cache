import os
import numpy as np
import h5py
import prody
from scipy.spatial import cKDTree
from scipy.ndimage import zoom
import subprocess
import shutil
from multiprocessing import Pool
import prody as pr


# ============================================================================
# RESOLUTION FIX: Fiziksel boyutu 161 Å olarak sabit tutuyoruz
# ============================================================================
BOX_SIZE_DEFAULT = 72
PHYSICAL_EXTENT = 161.0  # Ångström - 161x161x161 ile aynı fiziksel alan
APBS_SIZE_OFFSET = 1

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

def mol2_to_pdb(mol2_file, pdb_file, timeout=30):
    cmd = ["obabel", "-imol2", mol2_file, "-opdb", "-O", pdb_file]
    try:
        result = subprocess.run(cmd, capture_output=False, timeout=timeout, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.TimeoutExpired:
        #print(f"Timeout: OpenBabel conversion took too long for {mol2_file}")
        raise RuntimeError(f"obabel conversion timed out for {mol2_file}")
    if result.returncode != 0:
        #print("OpenBabel error:", result.stderr.decode())
        raise RuntimeError(f"obabel conversion failed for {mol2_file}")

def get_grid_params(protein_pdb, box_size, grid_step=1.0):
    """
    Grid parametrelerini hesapla.
    
    CHANGED: grid_step artık 161.0/box_size olarak otomatik hesaplanıyor
    Böylece fiziksel boyut 161 Å'de sabit kalıyor.
    """
    atoms = prody.parsePDB(protein_pdb).getCoords()
    center = atoms.mean(axis=0)
    grid_min = center - (box_size // 2) * grid_step
    grid_max = center + (box_size // 2) * grid_step
    grid_shape = (box_size, box_size, box_size)
    return grid_shape, grid_min, grid_max, center

def ligand_to_mask(ligand_pdb, grid_center, box_size, resolution):
    mol = prody.parsePDB(ligand_pdb)
    coords = mol.getCoords()
    grid_shape = (box_size, box_size, box_size)
    mask = np.zeros(grid_shape, dtype=np.uint8)
    grid_start = np.array(grid_center) - (box_size // 2) * resolution
    for atom in coords:
        idx = np.round((atom - grid_start) / resolution).astype(int)
        if np.all((idx >= 0) & (idx < box_size)):
            mask[tuple(idx)] = 1
    return mask

def prepare_protein_pocket(protein_pdb, ligand_pdb, pocket_pdb, selected_pdb):
    """
    Protein ve ligand PDB dosyalarini kullanarak ligand çevresindeki pocket'i seçip,
    temizlenmiş (h2o'suz) protein+ligand pocket pdb'lerini kaydeder.
    """
    # 1. Protein sadece protein atomlarını içeriyor olmalı
    structure = pr.parsePDB(protein_pdb)
    structure = structure.select('protein')
    pr.writePDB(protein_pdb, structure)

    # 2. Ligand PDB zaten hazır olmalı
    # 3. Liganda yakın zincirleri seç (exwithin 7A, tüm zincirler)
    protein = pr.parsePDB(protein_pdb)
    ligand = pr.parsePDB(ligand_pdb)
    lresname = ligand.getResnames()[0]
    complx = ligand + protein
    complx = complx.select(f'same chain as exwithin 7 of resname {lresname}')
    # Sadece protein atomlarını seç (ligand hariç)
    structure = complx.select(f'protein and not resname {lresname}')
    pr.writePDB(selected_pdb, structure)
    # ProDy header'ı kaldır
    with open(selected_pdb, 'r') as fin:
        data = fin.read().splitlines(True)
    with open(selected_pdb, 'w') as fout:
        fout.writelines(data[1:])
    # Pocket atomlarını seç (ligand çevresi, exwithin 4.5A)
    selected = pr.parsePDB(selected_pdb)
    complx = ligand + selected
    pocket = complx.select(f'same residue as exwithin 4.5 of resname {lresname}')
    pr.writePDB(pocket_pdb, pocket)
    return selected_pdb, pocket_pdb

def pdb2pqr_wrapper(selected_pdb, pqr_file, pdb2pqr_bin="pdb2pqr30"):
    cmd = [
        pdb2pqr_bin,
        "--with-ph=7.4",
        "--ff=AMBER",  # veya PARSE, iş akışına göre
        "--keep-chain",
        selected_pdb,
        pqr_file,
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    if proc.returncode != 0:
        print("pdb2pqr failed:", err.decode())
        raise RuntimeError(f"pdb2pqr failed: {err.decode()}")
    return pqr_file

def run_apbs(protein_pdb, protein_pqr, cache_dir, case_name, box_size, grid_step=1.0, timeout=60):
    """
    CHANGED: APBS her zaman 161x161x161 grid'de çalışıyor (spacing=1.0 Å).
    Sonra istenen box_size'a resize ediliyor.
    Bu sayede APBS dime kısıtlamaları (2^n*m+1 formatı) sorun çıkarmıyor.
    """
    apbs_case_dir = os.path.join(cache_dir, case_name)
    os.makedirs(apbs_case_dir, exist_ok=True)
    shutil.copy2(protein_pqr, apbs_case_dir)
    
    apbs_in_file = os.path.join(apbs_case_dir, "apbs.in")
    apbs_out_file = os.path.join(apbs_case_dir, "apbs_grid.dx")
    
    # CHANGED: APBS grid'i her zaman PHYSICAL_EXTENT ile 1.0 Å spacing'de hesapla
    apbs_box_size = int(PHYSICAL_EXTENT)  # 161
    apbs_grid_step = 1.0
    grid_shape, grid_min, grid_max, grid_center = get_grid_params(protein_pdb, apbs_box_size, apbs_grid_step)
    cglen = grid_max - grid_min
    dime = [apbs_box_size] * 3  # 161x161x161 - APBS için valid
    
    # FIXED: APBS için doğru grid spacing hesaplama
    # glen = fiziksel boyut, grid = voxel spacing
    grid_spacing = cglen[0] / (dime[0] - 1)  # Å/voxel
    
    with open(apbs_in_file, "w") as f:
        f.write(f"""read
    mol pqr {os.path.basename(protein_pqr)}
end
elec
    mg-manual
    mol 1
    dime {dime[0]} {dime[1]} {dime[2]}
    glen {cglen[0]:.3f} {cglen[1]:.3f} {cglen[2]:.3f}
    gcent {grid_center[0]:.3f} {grid_center[1]:.3f} {grid_center[2]:.3f}
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
        # DEBUG: APBS hatasını görmek için stdout/stderr'ı gösteriyoruz
        result = subprocess.run(["apbs", apbs_in_file], cwd=apbs_case_dir, check=True, timeout=timeout, 
                               capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"\n{'='*80}")
        print(f"APBS FAILED for {case_name}")
        print(f"Exit code: {e.returncode}")
        print(f"APBS input file: {apbs_in_file}")
        print(f"STDOUT:\n{e.stdout}")
        print(f"STDERR:\n{e.stderr}")
        print(f"{'='*80}\n")
        raise RuntimeError(f"apbs failed with exit code {e.returncode}: {e.stderr[:200]}")
    except subprocess.TimeoutExpired:
        print(f"Timeout: apbs took too long for {protein_pqr}")
        raise RuntimeError(f"apbs timed out for {protein_pqr}")
    apbs_grid = load_dx_grid(apbs_out_file)
    # CHANGED: APBS output'u (161x161x161) istenen box_size'a resize et
    apbs_grid_resized = resize_grid(apbs_grid, (box_size, box_size, box_size))
    
    # CHANGED: grid_center ve min/max'i son grid (box_size) için döndür
    final_grid_shape, final_grid_min, final_grid_max, _ = get_grid_params(protein_pdb, box_size, grid_step)
    return apbs_grid_resized, grid_center, final_grid_shape, final_grid_min, final_grid_max

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

def pdb_to_mask(pdb_file, grid_center, box_size, resolution):
    mol = prody.parsePDB(pdb_file)
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

def compute_distance_to_surface(protein_pdb, grid_center, box_size, resolution):
    mol = prody.parsePDB(protein_pdb)
    coords = mol.getCoords()
    grid = np.indices((box_size, box_size, box_size)).reshape(3, -1).T
    grid_start = np.array(grid_center) - (box_size // 2) * resolution
    grid_points = grid_start + grid * resolution
    tree = cKDTree(coords)
    dists, _ = tree.query(grid_points)
    dist_grid = dists.reshape((box_size, box_size, box_size)).astype(np.float32)
    return dist_grid

def process_pdbbind_case(case_dir, output_h5, box_size, resolution, output_dir=None, case_name=None, pdb2pqr_bin="pdb2pqr30"):
    protein_pdb = os.path.join(case_dir, f"{case_name}_protein.pdb")
    ligand_mol2 = os.path.join(case_dir, f"{case_name}_ligand.mol2")
    site_pdb = os.path.join(case_dir, f"{case_name}_pocket.pdb")    
    ligand_pdb = os.path.join(case_dir, f"{case_name}_ligand.pdb")
    selected_pdb = os.path.join(case_dir, f"{case_name}_selected.pdb")
    protein_pqr = os.path.join(case_dir, f"{case_name}_selected.pqr")
    site_calculated_pdb = os.path.join(case_dir, f"{case_name}_calculated_pocket.pdb")
    
    # 1. Ligand .mol2 -> .pdb çevir
    mol2_to_pdb(ligand_mol2, ligand_pdb)
    
    # 2. Liganda yakın protein pocket'ını hazırla (watersız ve zincir seçilmiş)
    selected_pdb, pocket_calculated_pdb = prepare_protein_pocket(protein_pdb, ligand_pdb, site_calculated_pdb, selected_pdb)

    # 3. pdb2pqr ile charge/radius ekle
    pdb2pqr_wrapper(selected_pdb, protein_pqr, pdb2pqr_bin)
    
      # 4. APBS gridini pocket üzerinden hesapla
    apbs_grid, grid_center, grid_shape, grid_min, grid_max = run_apbs(protein_pdb,
        protein_pqr, output_dir, case_name, box_size=box_size, grid_step=resolution
    )
  
    featurizer = Featurizer()
    ligand_mask = ligand_to_mask(ligand_pdb, grid_center, box_size, resolution)
    atomic_features = featurizer.featurize(protein_pdb, box_size, resolution, grid_center)
    hydro_grid = compute_hydrophobicity_grid(protein_pdb, grid_center, box_size, resolution)
    shape_mask = get_protein_shape_mask(protein_pdb, grid_center, box_size, resolution)
    dist2lig = compute_distance_to_ligand(ligand_pdb, grid_center, box_size, resolution)
    dist2surf = compute_distance_to_surface(protein_pdb, grid_center, box_size, resolution)
    label_calculated = pdb_to_mask(pocket_calculated_pdb, grid_center, box_size, resolution)
    label_in_dataset = pdb_to_mask(site_pdb, grid_center, box_size, resolution)

    
    with h5py.File(output_h5, "w") as h5f:
        save_atomic_features_separate(h5f, atomic_features, featurizer.FEATURE_NAMES)
        h5f.create_dataset("features/ligand", data=ligand_mask, compression="gzip")
        h5f.create_dataset("features/electrostatic_grid", data=apbs_grid, compression="gzip")
        h5f.create_dataset("features/shape", data=shape_mask, compression="gzip")
        h5f.create_dataset("features/dist_to_ligand", data=dist2lig, compression="gzip")
        h5f.create_dataset("features/hydrophobicity", data=hydro_grid, compression="gzip")
        h5f.create_dataset("features/dist_to_surface", data=dist2surf, compression="gzip")
        h5f.create_dataset("label/binding_site_calculated", data=label_calculated, compression="gzip")
        h5f.create_dataset("label/binding_site_in_dataset", data=label_in_dataset, compression="gzip")
        h5f.attrs["box_size"] = box_size
        h5f.attrs["resolution"] = resolution
        h5f.attrs["center"] = grid_center.tolist()
        # CHANGED: physical_extent bilgisini de kaydediyoruz
        h5f.attrs["physical_extent"] = PHYSICAL_EXTENT
        h5f.attrs["features"] = ",".join(
            [f"atomic_{x}" for x in featurizer.FEATURE_NAMES] +
            ["electrostatic_grid", "shape", "dist_to_ligand", "hydrophobicity", "dist_to_surface", "ligand"]
        )
        h5f.attrs["labels"] = ["binding_site_calculated", "binding_site_in_dataset"]
        h5f.attrs["electrostatic_grid_shape"] = [box_size, box_size, box_size]
        h5f.attrs["electrostatic_grid_min"] = grid_min.tolist()
        h5f.attrs["electrostatic_grid_max"] = grid_max.tolist()
    try:
        for f in [ligand_pdb]:
            if os.path.exists(f):
                os.remove(f)
        apbs_case_dir = os.path.join(output_dir, case_name)
        if os.path.exists(apbs_case_dir):
            shutil.rmtree(apbs_case_dir)
    except Exception as cleanup_err:
        print(f"Warning: For {case_dir} temporary files couldn't be deleted. Error message: {cleanup_err}")

def process_scpdb_case(case_dir, output_h5, box_size, resolution, output_dir=None, case_name=None, pdb2pqr_bin="pdb2pqr30"):
    protein_mol2 = os.path.join(case_dir, "protein.mol2")
    ligand_mol2 = os.path.join(case_dir, "ligand.mol2")
    site_mol2 = os.path.join(case_dir, "site.mol2")

    protein_pdb = os.path.join(case_dir, "protein.pdb")
    ligand_pdb = os.path.join(case_dir, "ligand.pdb")
    site_pdb = os.path.join(case_dir, "site.pdb")
    site_calculated_pdb = os.path.join(case_dir, f"{case_name}_site_calculated.pdb")
    selected_pdb = os.path.join(case_dir, "protein_selected.pdb")
    protein_pqr = os.path.join(case_dir, "protein_selected.pqr")

    # 1. Mol2 dosyalarını pdb'ye çevir
    mol2_to_pdb(protein_mol2, protein_pdb)
    mol2_to_pdb(ligand_mol2, ligand_pdb)
    mol2_to_pdb(site_mol2, site_pdb)
    
      # 2. Protein pocket'ını hazırla (genellikle site pocket ile aynı)
    selected_pdb, pocket_calculated_pdb = prepare_protein_pocket(protein_pdb, ligand_pdb, site_calculated_pdb, selected_pdb)

    # 3. pdb2pqr ile charge/radius ekle
    pdb2pqr_wrapper(selected_pdb, protein_pqr, pdb2pqr_bin)


    # 4. APBS gridini protein pocket üzerinden hesapla
    apbs_grid, grid_center, grid_shape, grid_min, grid_max = run_apbs(protein_pdb,
        protein_pqr, output_dir, case_name, box_size=box_size, grid_step=resolution
    )

    featurizer = Featurizer()
    ligand_mask = ligand_to_mask(ligand_pdb, grid_center, box_size, resolution)
    atomic_features = featurizer.featurize(protein_pdb, box_size, resolution, grid_center)
    hydro_grid = compute_hydrophobicity_grid(protein_pdb, grid_center, box_size, resolution)
    shape_mask = get_protein_shape_mask(protein_pdb, grid_center, box_size, resolution)
    dist2lig = compute_distance_to_ligand(ligand_pdb, grid_center, box_size, resolution)
    label = pdb_to_mask(pocket_calculated_pdb, grid_center, box_size, resolution)
    dist2surf = compute_distance_to_surface(protein_pdb, grid_center, box_size, resolution)
    label_in_dataset = pdb_to_mask(site_pdb, grid_center, box_size, resolution)

    with h5py.File(output_h5, "w") as h5f:
        save_atomic_features_separate(h5f, atomic_features, featurizer.FEATURE_NAMES)
        h5f.create_dataset("features/ligand", data=ligand_mask, compression="gzip")    
        h5f.create_dataset("features/electrostatic_grid", data=apbs_grid, compression="gzip")
        h5f.create_dataset("features/shape", data=shape_mask, compression="gzip")
        h5f.create_dataset("features/dist_to_ligand", data=dist2lig, compression="gzip")
        h5f.create_dataset("features/hydrophobicity", data=hydro_grid, compression="gzip")
        h5f.create_dataset("features/dist_to_surface", data=dist2surf, compression="gzip")
        h5f.create_dataset("label/binding_site_calculated", data=label, compression="gzip")
        h5f.create_dataset("label/binding_site_in_dataset", data=label_in_dataset, compression="gzip")
        h5f.attrs["box_size"] = box_size
        h5f.attrs["resolution"] = resolution
        h5f.attrs["center"] = grid_center.tolist()
        # CHANGED: physical_extent bilgisini de kaydediyoruz
        h5f.attrs["physical_extent"] = PHYSICAL_EXTENT
        h5f.attrs["features"] = ",".join(
            [f"atomic_{x}" for x in featurizer.FEATURE_NAMES] +
            ["electrostatic_grid", "shape", "", "hydrophobicity", "dist_to_surface", "ligand"]
        )
        h5f.attrs["labels"] = ["binding_site_calculated", "binding_site_in_dataset"]
        h5f.attrs["electrostatic_grid_shape"] = [box_size, box_size, box_size]
        h5f.attrs["electrostatic_grid_min"] = grid_min.tolist()
        h5f.attrs["electrostatic_grid_max"] = grid_max.tolist()
    
    try:
        for f in [protein_pdb, ligand_pdb, site_pdb]:
            if os.path.exists(f):
                os.remove(f)
        apbs_case_dir = os.path.join(output_dir, case_name)
        if os.path.exists(apbs_case_dir):
            shutil.rmtree(apbs_case_dir)
    except Exception as cleanup_err:
        print(f"Warning: For {case_dir} temporary files couldn't be deleted. Error message: {cleanup_err}")


def process_case(args):
    case_dir, out_h5, box_size, resolution, output_dir, case_name, dataset_type, failed_cases_path = args
    output_protein_dir = os.path.join(output_dir, case_name)
    if os.path.exists(out_h5):
        print(f"[SKIP] Already exists: {out_h5}")
        return
    try:
        if dataset_type == "pdbbind":
            process_pdbbind_case(
                case_dir, out_h5, box_size=box_size, resolution=resolution, output_dir=output_dir, case_name=case_name
            )
        elif dataset_type == "scpdb":
            process_scpdb_case(
                case_dir, out_h5, box_size=box_size, resolution=resolution, output_dir=output_dir, case_name=case_name
            )
        else:
            raise ValueError("Unknown dataset type: " + str(dataset_type))
    except Exception as e:
        print(f"[FAILED] {case_name}: {e}")
        # Hatalı klasörü dosyaya yaz
        with open(failed_cases_path, "a") as f:
            f.write(case_name + "\n")            
        # Delete problematic input folder
        try:
            print(f"Deleting failed case folder: {output_protein_dir}")
            # shutil.rmtree(output_protein_dir)
        except Exception as cleanup_err:
            with open(failed_cases_path, "a") as f:
                f.write(f"Error deleting folder {output_protein_dir}: {cleanup_err}\n")
        # Delete output .h5 if exists (partial/garbage file)
        try:
            if os.path.exists(out_h5):
                print(f"Deleting failed output: {out_h5}")
                os.remove(out_h5)
        except Exception as cleanup_err:
            with open(failed_cases_path, "a") as f:
                f.write(f"Error deleting output {out_h5}: {cleanup_err}\n")

def get_dataset_paths(dataset="pdbbind", box_size=72):
    """
    CHANGED: Output directory'ye "_resolution_fix" suffix eklendi
    Böylece eski cache'lerle karışmıyor.
    """
    if dataset == "pdbbind":
        pdbbind_root = "/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set/"
        output_dir = f"/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set_minimal_cache_only_fits_60_resolution_fix/box{box_size}/"
    elif dataset == "scpdb":
        pdbbind_root = "/Users/tevfik/Sandbox/github/PHD/data/scPDB/"
        output_dir = f"/Users/tevfik/Sandbox/github/PHD/data/scPDB_cache_only_fits_60_resolution_fix/box{box_size}/"
    else:
        raise ValueError("dataset must be 'pdbbind' or 'scpdb'")
    return pdbbind_root, output_dir

if __name__ == "__main__":
    import argparse
    import tqdm

    parser = argparse.ArgumentParser(description="Generate cache with RESOLUTION FIX - maintains 161 Å physical extent")
    parser.add_argument("--dataset", type=str, choices=["pdbbind", "scpdb"], default="pdbbind")
    parser.add_argument("--box_size", type=int, default=BOX_SIZE_DEFAULT,
                        help="Grid size in voxels (default: 72). Physical extent stays at 161 Å")
    # CHANGED: resolution parametresi artık PHYSICAL_EXTENT/box_size olarak otomatik hesaplanıyor
    # Kullanıcı manuel olarak override edebilir ama önerilmiyor
    parser.add_argument("--resolution", type=float, default=None,
                        help="Voxel size in Å (default: auto-calculated as 161.0/box_size)")
    parser.add_argument("--nproc", type=int, default=6)
    parser.add_argument("--max_size", type=int, default=60)
    args = parser.parse_args()

    # CHANGED: Resolution'ı otomatik hesapla
    if args.resolution is None:
        args.resolution = PHYSICAL_EXTENT / args.box_size
        print(f"Auto-calculated resolution: {args.resolution:.4f} Å/voxel")
        print(f"Physical extent: {PHYSICAL_EXTENT} Å × {PHYSICAL_EXTENT} Å × {PHYSICAL_EXTENT} Å")
        print(f"Grid size: {args.box_size} × {args.box_size} × {args.box_size} voxels")
    else:
        print(f"WARNING: Using manual resolution {args.resolution} Å/voxel")
        print(f"Physical extent will be: {args.resolution * args.box_size:.1f} Å (not 161 Å!)")

    dataset_root, output_dir = get_dataset_paths(args.dataset, args.box_size)
    os.makedirs(output_dir, exist_ok=True)
    failed_cases_path = os.path.join(output_dir, "failed_cases.txt")
    # Eski log'u sil
    if os.path.exists(failed_cases_path):
        os.remove(failed_cases_path)

    #all_cases = sorted([f for f in os.listdir(pdbbind_root) if os.path.isdir(os.path.join(pdbbind_root, f))])
    dataset_fits_file = f"{args.dataset}_{args.max_size}_fits.txt"
    with open(dataset_fits_file, 'r') as f:
        
        all_cases = [line.strip() for line in f if line.strip()]
    print(f"Processing {len(all_cases)} cases in {args.dataset} with box size {args.box_size} and {args.nproc} processes.")

    tasks = []
    for case in all_cases:
        case_dir = os.path.join(dataset_root, case)
        out_h5 = os.path.join(output_dir, f"{case}.h5")
        tasks.append((case_dir, out_h5, args.box_size, args.resolution, output_dir, case, args.dataset, failed_cases_path))

    with Pool(processes=args.nproc) as pool:
        list(tqdm.tqdm(pool.imap_unordered(process_case, tasks), total=len(tasks)))

    print("ALL DONE.")
    print(f"Failed cases written to: {failed_cases_path}")
