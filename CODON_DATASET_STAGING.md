# Codon dataset staging and Work11 cache generation

This workflow stages the raw datasets under the shared NFS root and writes the
generated H5 cache files under the same root.

## Root

```text
/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS
```

Expected layout after staging:

```text
DEEP_APBS_DATASETS/
  archives/
  staging/
  datasets/
    scPDB/
    pdbbind/
      refined-set -> ../pdbbind2020/refined-set
    pdbbind2020/
      refined-set/
    external_benchmarks/
      puresnet/
      p2rank-datasets/
  cache/
    work11_cache_gridfix_v1/
```

## Stage datasets on datamover

```bash
cd /homes/tevfik/PHD/generate_cache
git pull

DATA_ROOT=/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS \
sbatch scripts/run_codon_datamover_prepare_datasets.sbatch
```

Monitor:

```bash
squeue -u "$USER"
tail -f /homes/tevfik/PHD/slurm_logs/deep-apbs-stage-data-<jobid>.out
```

## Transfer local datasets to Codon NFS

If the local workstation already has PDBBind and external benchmark files, use
the rsync helper from the local machine. The remote rsync process is launched
inside the Codon `datamover` partition, so writes to `/nfs/production/...` do
not happen on the login node.

Dry run:

```bash
cd /Users/tevfik/Sandbox/Tevfik/Projects/phd_examples/generate_cache

DRY_RUN=1 \
scripts/rsync_local_datasets_to_codon_datamover.sh
```

Actual transfer:

```bash
cd /Users/tevfik/Sandbox/Tevfik/Projects/phd_examples/generate_cache

scripts/rsync_local_datasets_to_codon_datamover.sh
```

Default local inputs:

```text
/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set
/Users/tevfik/Sandbox/github/PHD/data/external_benchmarks
```

Default remote outputs:

```text
/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS/datasets/pdbbind/refined-set
/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS/datasets/pdbbind2020/refined-set
/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS/datasets/external_benchmarks
```

## PDBBind note

PDBBind usually requires licensed/manual access. The staging script therefore
expects the official archive to be present locally, unless an explicit
`PDBBIND2020_REFINED_URL` environment variable is provided.

Place the official archive in:

```text
/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS/archives/
```

Use one of these names:

```text
pdbbind_v2020_refined.tar.gz
PDBbind_v2020_refined.tar.gz
pdbbind_refined.tar.gz
```

Then rerun only the PDBBind staging step:

```bash
DATA_ROOT=/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS \
SKIP_SCPDB=1 \
SKIP_PURESNET=1 \
SKIP_P2RANK=1 \
FORCE=1 \
sbatch scripts/run_codon_datamover_prepare_datasets.sbatch
```

If a direct refined-set URL is available in your PDBBind account/session:

```bash
PDBBIND2020_REFINED_URL=https://.../pdbbind_v2020_refined.tar.gz \
DATA_ROOT=/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS \
SKIP_SCPDB=1 \
SKIP_PURESNET=1 \
SKIP_P2RANK=1 \
sbatch scripts/run_codon_datamover_prepare_datasets.sbatch
```

If a separate non-2020 PDBBind refined-set URL is available:

```bash
PDBBIND_REFINED_URL=https://.../pdbbind_refined.tar.gz \
DATA_ROOT=/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS \
sbatch scripts/run_codon_datamover_prepare_datasets.sbatch
```

## Smoke cache run

Run a small check before launching the full cache matrix:

```bash
DATA_ROOT=/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS \
GRID_SPECS="36:70" \
LIMIT=20 \
NPROC=8 \
sbatch scripts/run_work11_slurm_cache_matrix_nfs.sbatch
```

Monitor:

```bash
tail -f /homes/tevfik/PHD/slurm_logs/work11-cache-nfs-<jobid>.out
```

## Full cache matrix

```bash
DATA_ROOT=/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS \
GRID_SPECS="36:70 72:160 161:160" \
LIMIT=all \
NPROC=16 \
sbatch scripts/run_work11_slurm_cache_matrix_nfs.sbatch
```

Outputs:

```text
/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS/cache/work11_cache_gridfix_v1/scpdb/label_cavity6/box36_span70/
/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS/cache/work11_cache_gridfix_v1/scpdb/label_cavity6/box72_span160/
/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS/cache/work11_cache_gridfix_v1/scpdb/label_cavity6/box161_span160/
/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS/cache/work11_cache_gridfix_v1/pdbbind/refined-set/box36_span70/
/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS/cache/work11_cache_gridfix_v1/pdbbind/refined-set/box72_span160/
/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS/cache/work11_cache_gridfix_v1/pdbbind/refined-set/box161_span160/
```
