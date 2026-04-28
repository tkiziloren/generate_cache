import os
import subprocess
import prody
from multiprocessing import Pool
from tqdm import tqdm

# --- Ayarlar ---
BOX_SIZE      = 128
PDBBIND_ROOT  = "/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set"
SCPDB_ROOT    = "/Users/tevfik/Sandbox/github/PHD/data/scPDB"
MAX_PROCESSES = 16

def mol2_to_pdb_obabel(mol2_file, pdb_file, timeout=30):
    cmd = ["obabel", "-imol2", mol2_file, "-opdb", "-O", pdb_file]
    subprocess.run(cmd, check=True, timeout=timeout,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def get_dims_from_pdb(pdb_file):
    coords = prody.parsePDB(pdb_file).getCoords()
    return coords.max(axis=0) - coords.min(axis=0)

def pdbbind_worker(case):
    pdb_file = os.path.join(PDBBIND_ROOT, case, f"{case}_protein.pdb")
    dims     = get_dims_from_pdb(pdb_file)
    max_dim  = dims.max()
    req = 0 if max_dim <= BOX_SIZE else int(max_dim) + (0 if max_dim.is_integer() else 1)
    return case, req

def scpdb_worker(case):
    case_dir = os.path.join(SCPDB_ROOT, case)
    mol2     = os.path.join(case_dir, "protein.mol2")
    tmp_pdb  = os.path.join(case_dir, "protein_tmp.pdb")
    mol2_to_pdb_obabel(mol2, tmp_pdb)
    dims    = get_dims_from_pdb(tmp_pdb)
    os.remove(tmp_pdb)
    max_dim = dims.max()
    req = 0 if max_dim <= BOX_SIZE else int(max_dim) + (0 if max_dim.is_integer() else 1)
    return case, req

def run_check(root, cases, worker, desc):
    valid, fits, unfits = [], [], []
    for c in cases:
        if c.startswith("."):
            continue
        path = os.path.join(root, c, f"{c}_protein.pdb") if worker is pdbbind_worker else os.path.join(root, c, "protein.mol2")
        if os.path.isfile(path):
            valid.append(c)
        else:
            tqdm.write(f"[SKIP {desc}] {c}: {os.path.basename(path)} yok")
    with Pool(processes=MAX_PROCESSES) as pool:
        for case, req in tqdm(pool.imap_unordered(worker, valid),
                              total=len(valid), desc=desc, unit="case"):
            if req == 0:
                fits.append(case)
            else:
                unfits.append((case, req))
    # sort unfits by required grid descending
    unfits.sort(key=lambda x: x[1], reverse=True)
    return fits, unfits

if __name__ == "__main__":
    # PDBbind kontrolü
    pb_cases = sorted(os.listdir(PDBBIND_ROOT))
    fits_pb, unfits_pb = run_check(PDBBIND_ROOT, pb_cases, pdbbind_worker, "PDBbind")

    # scPDB kontrolü
    sc_cases = sorted(os.listdir(SCPDB_ROOT))
    fits_sc, unfits_sc = run_check(SCPDB_ROOT, sc_cases, scpdb_worker, "scPDB")

    # Dosyalara yazma
    box_tag = str(int(BOX_SIZE))
    with open(f"pdbbind_{box_tag}_fits.txt", "w") as f:
        for case in fits_pb:
            f.write(f"{case}\n")
    with open(f"pdbbind_{box_tag}_unfits.txt", "w") as f:
        for case, req in unfits_pb:
            f.write(f"{case} requires {req}Å\n")

    with open(f"scpdb_{box_tag}_fits.txt", "w") as f:
        for case in fits_sc:
            f.write(f"{case}\n")
    with open(f"scpdb_{box_tag}_unfits.txt", "w") as f:
        for case, req in unfits_sc:
            f.write(f"{case} requires {req}Å\n")

    # Konsola özet
    print(f"\n--- PDBbind refined-set ({len(fits_pb)+len(unfits_pb)} cases) ---")
    print(f"Fits in {BOX_SIZE}Å: {len(fits_pb)}")
    print(f"Does NOT fit : {len(unfits_pb)}")
    print("  Unfits (largest→smallest):")
    for case, req in unfits_pb:
        print(f"    {case}: requires {req}Å")

    print(f"\n--- scPDB ({len(fits_sc)+len(unfits_sc)} cases) ---")
    print(f"Fits in {BOX_SIZE}Å: {len(fits_sc)}")
    print(f"Does NOT fit : {len(unfits_sc)}")
    print("  Unfits (largest→smallest):")
    for case, req in unfits_sc:
        print(f"    {case}: requires {req}Å")