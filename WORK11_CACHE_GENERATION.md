# Work11 - 36/72/161 Gridfix Cache Generation

Bu work package tez taslagindan bagimsizdir. Amac, scPDB ve PDBBind icin ayni H5 semasinda 36, 72 ve 161 grid cache uretimini temiz ve tekrar calistirilabilir hale getirmektir.

## Neden Work11?

Mevcut 36 grid scPDB uretim loglarinda hata grid boyutundan degil, PDB2PQR/APBS hazirlik adimindan geliyordu. Tipik hata:

```text
CRITICAL:Unable to debump biomolecule. Biomolecular structure is incomplete
Valist_readPQR: Error parsing atom
```

Bu hatalar cogunlukla scPDB `protein.mol2` dosyalarindan PDB'ye gecen hidrojen atomlarinin veya NAD gibi kofaktorlerin PDB2PQR tarafinda yeniden protonlama/debump asamasini bozmasindan kaynaklaniyor. Work11 icin gridfix cache kodu artik PDB2PQR'a ve protein feature uretimine sadece standart amino-acid heavy atom setini verir:

```text
protein_atom_policy = standard_amino_acid_heavy_atoms_only
```

Bu degisiklik hem scPDB hem de PDBBind gridfix cache akisi icin gecerlidir.

## H5 Dosya Icerigi

Her protein/case icin tek bir `.h5` dosyasi uretilir.

Genel yapi:

```text
case.h5
  attrs/
  features/
  label/
```

### Ortak Attributes

Yeni gridfix H5 dosyalarinda beklenen ana metadata:

```text
schema_version
grid_convention
box_size
resolution
center
grid_origin
grid_max
physical_span_angstrom
apbs_grid_points
apbs_span_angstrom
electrostatic_grid_shape
electrostatic_grid_min
electrostatic_grid_max
features
labels
protein_atom_policy
```

scPDB icin ek metadata:

```text
dataset = scpdb
case_name
dataset_label_source
include_all_labels
all_label_sources
requested_target_span_angstrom
```

### Features

Tum feature tensorleri ayni shape'tedir:

```text
(box_size, box_size, box_size)
```

Su an yazdigimiz featurelar:

```text
features/atomic_B
features/atomic_C
features/atomic_N
features/atomic_O
features/atomic_P
features/atomic_S
features/atomic_Se
features/atomic_halogen
features/atomic_metal
features/atomic_hyb
features/atomic_heavydegree
features/atomic_heterodegree
features/atomic_partialcharge
features/atomic_molcode
features/atomic_hydrophobic
features/atomic_aromatic
features/atomic_acceptor
features/atomic_donor
features/atomic_ring
features/ligand
features/electrostatic_grid
features/shape
features/dist_to_ligand
features/hydrophobicity
features/dist_to_surface
```

Notlar:

- `electrostatic_grid` APBS potansiyel grididir. APBS kaynak grid 161 nokta ve 160 A span ile uretilir, sonra hedef grid boyutuna resample edilir.
- `shape` standart amino-acid heavy atom maskesidir.
- `ligand` ligand atom maskesidir.
- `dist_to_ligand` deployment icin leakage riskli bir feature'dir; cache icinde analiz/label kontrolu icin durabilir ama final model input kombinasyonlarina dahil edilmemelidir.
- `dist_to_surface` ve `hydrophobicity` leakage degildir; yine de feature ablation ile etkileri ayri olculmelidir.

### Legacy + V2 Feature Extension

Work11 sadece eski cache'i yeniden uretmeyecek. Eski feature isimleri aynen korunacak, ayni H5 dosyasina iyilestirilmis APBS ve atomic feature kanallari da yeni isimlerle eklenecek. Bunun nedeni eski deneylerin tekrar uretilebilir kalmasi ve yeni representation etkisinin kontrollu ablation ile olculebilmesidir.

Legacy feature isimleri degistirilmeyecek:

```text
features/electrostatic_grid
features/shape
features/ligand
features/dist_to_ligand
features/hydrophobicity
features/dist_to_surface
features/atomic_B
features/atomic_C
features/atomic_N
features/atomic_O
features/atomic_P
features/atomic_S
features/atomic_Se
features/atomic_halogen
features/atomic_metal
features/atomic_hyb
features/atomic_heavydegree
features/atomic_heterodegree
features/atomic_partialcharge
features/atomic_molcode
features/atomic_hydrophobic
features/atomic_aromatic
features/atomic_acceptor
features/atomic_donor
features/atomic_ring
```

V2 APBS kanallari:

```text
features/electrostatic_grid_v2_full_protein_raw
features/electrostatic_grid_v2_full_protein_clip20_minmax
features/electrostatic_grid_v2_full_protein_signed
features/electrostatic_grid_v2_full_protein_pos
features/electrostatic_grid_v2_full_protein_neg
features/electrostatic_grid_v2_selected_chains_raw
```

V2 APBS atom policy:

```text
protein_original
-> protein_clean_full          # ligand-free, cleaned full protein context
-> protein_apbs_full           # default APBS v2 input
-> protein_apbs_selected       # legacy/current comparison APBS input
-> protein_pocket_label        # calculated binding-site label only
```

APBS v2 default'u `protein_apbs_full` olmalidir. `selected_chains` sadece legacy davranisi acik isimle tekrar olcmek icin tutulur.

V2 MOL2 atomic kanallari:

```text
features/atomic_mol2_B
features/atomic_mol2_C
features/atomic_mol2_N
features/atomic_mol2_O
features/atomic_mol2_P
features/atomic_mol2_S
features/atomic_mol2_Se
features/atomic_mol2_halogen
features/atomic_mol2_metal
features/atomic_mol2_hyb
features/atomic_mol2_heavydegree
features/atomic_mol2_heterodegree
features/atomic_mol2_partialcharge
features/atomic_mol2_molcode
features/atomic_mol2_hydrophobic
features/atomic_mol2_aromatic
features/atomic_mol2_acceptor
features/atomic_mol2_donor
features/atomic_mol2_ring
```

V2 geometry/context kanallari:

```text
features/shape_v2_full_protein
features/shape_v2_selected_chains
features/dist_to_surface_v2_full_protein
features/hydrophobicity_v2_full_protein
features/metal_mask_v2
features/modified_residue_mask_v2
features/cofactor_mask_v2
```

Leakage riski nedeniyle su kanallar final de novo prediction input'una girmemelidir:

```text
features/ligand
features/dist_to_ligand
features/dist_to_ligand_v2
features/dist_to_label_v2
```

### V2 Audit Metadata

V2 cache uretiminde H5 attrs ve manifest CSV icinde atom kaybini takip edecek alanlar bulunmalidir:

```text
protein_atom_policy_legacy
protein_atom_policy_v2
apbs_source_legacy
apbs_source_v2_full_protein
apbs_source_v2_selected_chains
atomic_source_legacy
atomic_source_v2_mol2
original_protein_atoms
clean_full_protein_atoms
selected_chain_atoms
apbs_full_protein_atoms
apbs_selected_chain_atoms
atomic_mol2_atoms
dropped_hydrogen_atoms
dropped_ligand_atoms
dropped_unknown_hetatm_atoms
included_chains
excluded_chains
apbs_pqr_atom_count
pdb2pqr_status
apbs_status
```

Bu alanlar olmadan v2 feature'larin eski feature'lardan neden farkli davrandigini yorumlamak zor olur.

### Labels

PDBBind:

```text
label/binding_site_calculated
label/binding_site_in_dataset
```

scPDB:

```text
label/binding_site_calculated
label/binding_site_in_dataset
```

`--include-all-labels` kullanildiginda scPDB icin ek olarak:

```text
label/binding_site_site
label/binding_site_cavity6
label/binding_site_cavityALL
```

scPDB'de `binding_site_in_dataset`, `--dataset-label-source` ile secilen label'in trainer uyumlu alias'idir. Work11 default'u:

```text
--dataset-label-source cavity6
--include-all-labels
```

## Grid Boyutlari

Work11 icin 36 gridde eski Kalasanty-benzeri deney hattini yeniden uretmek amaciyla 70 A span defaulttur. 72 ve 161 gridlerde APBS kaynak gridine uyumlu 160 A span kullanilir.

| box_size | target_span | voxel spacing |
|---:|---:|---:|
| 36 | 70 A | 70 / 35 = 2.0000 A |
| 72 | 160 A | 160 / 71 = 2.2535 A |
| 161 | 160 A | 160 / 160 = 1.0000 A |

APBS yine 161 x 161 x 161 kaynak grid ve 160 A span ile hesaplanir; hedef grid hangi span/resolution istiyorsa oraya resample edilir. 36 grid Kalasanty/PUResNet tarzina daha yakin hizli deney grididir; 72 ve 161 daha yuksek cozumlu ileri deneyler icindir.

## Local 36 Cache

Once lokalde 36 grid icin smoke/full uretim yapilacak.

Script:

```text
/Users/tevfik/Sandbox/Tevfik/Projects/phd_examples/generate_cache/scripts/run_work11_local_box36_cache.sh
```

Guvenli smoke test:

```bash
cd /Users/tevfik/Sandbox/Tevfik/Projects/phd_examples/generate_cache

LIMIT=20 \
NPROC=4 \
BOX_SIZE=36 \
TARGET_SPAN=70 \
LOCAL_OUTPUT_ROOT=/Users/tevfik/Sandbox/github/PHD/data/work11_local_cache_gridfix_v1 \
scripts/run_work11_local_box36_cache.sh
```

Local full 36 cache:

```bash
cd /Users/tevfik/Sandbox/Tevfik/Projects/phd_examples/generate_cache

LIMIT=all \
NPROC=8 \
BOX_SIZE=36 \
TARGET_SPAN=70 \
LOCAL_OUTPUT_ROOT=/Users/tevfik/Sandbox/github/PHD/data/work11_local_cache_gridfix_v1 \
scripts/run_work11_local_box36_cache.sh
```

Local ciktilar:

```text
/Users/tevfik/Sandbox/github/PHD/data/work11_local_cache_gridfix_v1/scpdb/label_cavity6/box36_span70/
/Users/tevfik/Sandbox/github/PHD/data/work11_local_cache_gridfix_v1/pdbbind/refined-set/box36_span70/
```

Kontrol:

```bash
tail -f /Users/tevfik/Sandbox/github/PHD/data/work11_local_cache_gridfix_v1/scpdb/label_cavity6/box36_span70/generation.log
tail -f /Users/tevfik/Sandbox/github/PHD/data/work11_local_cache_gridfix_v1/pdbbind/refined-set/box36_span70/generation.log
```

Basarisiz scPDB case'leri:

```text
scpdb/.../manifest.csv
```

Basarisiz PDBBind case'leri:

```text
pdbbind/.../failed_cases.txt
```

## Codon SLURM 36/72/161 Cache

Script:

```text
/Users/tevfik/Sandbox/Tevfik/Projects/phd_examples/generate_cache/scripts/run_work11_slurm_cache_matrix.sbatch
```

Codon komutu:

```bash
cd /homes/tevfik/PHD/phd_examples/generate_cache

GRID_SPECS="36:70 72:120 161:160" \
LIMIT=all \
NPROC=16 \
OUTPUT_ROOT=/hps/nobackup/arl/chembl/tevfik/deep-apbs-data/work11_cache_gridfix_v1 \
sbatch scripts/run_work11_slurm_cache_matrix.sbatch
```

Codon ciktilar:

```text
/hps/nobackup/arl/chembl/tevfik/deep-apbs-data/work11_cache_gridfix_v1/scpdb/label_cavity6/box36_span70/
/hps/nobackup/arl/chembl/tevfik/deep-apbs-data/work11_cache_gridfix_v1/scpdb/label_cavity6/box72_span120/
/hps/nobackup/arl/chembl/tevfik/deep-apbs-data/work11_cache_gridfix_v1/scpdb/label_cavity6/box161_span160/

/hps/nobackup/arl/chembl/tevfik/deep-apbs-data/work11_cache_gridfix_v1/pdbbind/refined-set/box36_span70/
/hps/nobackup/arl/chembl/tevfik/deep-apbs-data/work11_cache_gridfix_v1/pdbbind/refined-set/box72_span120/
/hps/nobackup/arl/chembl/tevfik/deep-apbs-data/work11_cache_gridfix_v1/pdbbind/refined-set/box161_span160/
```

SLURM loglari:

```text
/homes/tevfik/PHD/slurm_logs/work11-cache-<jobid>.out
/homes/tevfik/PHD/slurm_logs/work11-cache-<jobid>.err
```

## Uretim Sonrasi Kontrol

Her klasor icin kontrol edilecekler:

```bash
find <cache_dir> -maxdepth 1 -name "*.h5" | wc -l
tail -50 <cache_dir>/generation.log
```

scPDB:

```bash
python - <<'PY'
import pandas as pd
m = pd.read_csv("<scpdb_cache_dir>/manifest.csv")
print(m["status"].value_counts(dropna=False))
print(m[m["status"] == "failed"][["case", "error"]].head(30).to_string(index=False))
PY
```

PDBBind:

```bash
cat <pdbbind_cache_dir>/failed_cases.txt
```

## Beklenen Sonuc

1. Local 36 smoke/full cache ile PDB2PQR/APBS hatasinin azaldigi dogrulanacak.
2. Codon'da scPDB ve PDBBind icin 36, 72 ve 161 cache setleri uretilecek.
3. Her grid boyutu ayni H5 semasini kullanacagi icin training tarafinda sadece `h5_directory`, `box_size` ve feature set degisecek.
4. Basarisiz kalan case'ler manifest/failed_cases uzerinden kategorize edilip gerekiyorsa ikinci temizlik adimi planlanacak.

## Yapilan Smoke Test

2026-05-11 tarihinde iki kisa smoke test calisti:

```text
scPDB 1a4z_4, box36_span70 -> OK
PDBBind 1a1e, box36_span70 -> OK
```

`1a4z_4`, eski `box36_span70` logunda PDB2PQR hatasi veren case'lerden biriydi. Yeni standart amino-acid heavy atom filtresiyle H5 uretimi basarili oldu.
