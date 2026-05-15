# scPDB Gridfix Cache

This workflow generates one HDF5 file per scPDB case. It keeps the same cache schema used by the configurable 3D U-Net training code:

- `features/atomic_*`
- `features/ligand`
- `features/electrostatic_grid`
- `features/shape`
- `features/dist_to_ligand`
- `features/hydrophobicity`
- `features/dist_to_surface`
- `label/binding_site_calculated`
- `label/binding_site_in_dataset`

The scPDB root used by default is:

```bash
/Users/tevfik/Sandbox/github/PHD/data/scPDB
```

Each input case is expected to contain:

```text
protein.mol2
ligand.mol2
site.mol2
cavity6.mol2
cavityALL.mol2
```

## Label Sources

`binding_site_calculated` is generated from the ligand using the current protein-pocket selection logic.

`binding_site_in_dataset` is generated from one selected scPDB file:

- `--dataset-label-source site`
- `--dataset-label-source cavity6`
- `--dataset-label-source cavityALL`

Use `site` to match the current older cache behavior. Use `cavity6` or `cavityALL` when testing cavity-point labels closer to a scPDB/VolSite-style setup.

When `--include-all-labels` is set, every H5 also contains:

- `label/binding_site_site`
- `label/binding_site_cavity6`
- `label/binding_site_cavityALL`

`label/binding_site_in_dataset` remains as the trainer-compatible alias for the selected `--dataset-label-source`.

## Smoke Run

Generate a very small test cache:

```bash
python3 generate_cache_scpdb_gridfix.py \
  --dataset-label-source cavity6 \
  --include-all-labels \
  --box-size 72 \
  --target-span 160 \
  --limit 3 \
  --nproc 1
```

The output defaults to:

```bash
/Users/tevfik/Sandbox/github/PHD/data/scPDB_cache_gridfix_v1/label_cavity6/box72_span160
```

## Kalasanty-Style 36 Grid Test

For a 36-grid experiment, keep the physical span explicit. A span of `160 A` is directly comparable with the current gridfix APBS setup but coarse. A smaller span such as `70 A` is a separate local-grid experiment and must be validated for coverage before training.

```bash
python3 generate_cache_scpdb_gridfix.py \
  --dataset-label-source cavity6 \
  --include-all-labels \
  --box-size 36 \
  --target-span 160 \
  --limit 10 \
  --nproc 1
```

Alternative local-grid test:

```bash
python3 generate_cache_scpdb_gridfix.py \
  --dataset-label-source cavity6 \
  --include-all-labels \
  --box-size 36 \
  --target-span 70 \
  --limit 10 \
  --nproc 1
```

Do not compare `box36_span160` and `box36_span70` as the same experiment. They have different physical resolution and different coverage.

## Full Run

Generate all scPDB cases with the selected label source:

```bash
python3 generate_cache_scpdb_gridfix.py \
  --dataset-label-source cavity6 \
  --include-all-labels \
  --box-size 72 \
  --target-span 160 \
  --nproc 6
```

Each output directory contains `manifest.csv` with one row per case:

- status
- output H5 path
- resolution
- label source
- protein/ligand/label atom coverage inside the grid
- ligand voxel count
- calculated label voxel count
- dataset label voxel count
- per-source label coverage for `site`, `cavity6`, and `cavityALL` when `--include-all-labels` is used
- protein shape voxel count
- error message for failed cases

Use the manifest before training to check whether labels are non-empty and whether the chosen physical span covers the expected pocket.

The same output directory also contains `generation.log`. Use it for live monitoring and post-run analysis:

```bash
tail -f /Users/tevfik/Sandbox/github/PHD/data/scPDB_cache_gridfix_v1/label_cavity6/box36_span70/generation.log
```

## Kalasanty Folds

Kalasanty's official exclusion and fold files are stored under:

```bash
data/kalasanty
```

After cache generation finishes, build trainer-ready split lists:

```bash
PATH="$PWD/.conda/bin:$PATH" .conda/bin/python build_kalasanty_scpdb_splits.py \
  --h5-dir /Users/tevfik/Sandbox/github/PHD/data/scPDB_cache_gridfix_v1/label_cavity6/box36_span70 \
  --entry-policy all
```

This writes:

```bash
/Users/tevfik/Sandbox/github/PHD/data/scPDB_cache_gridfix_v1/label_cavity6/box36_span70/splits_kalasanty
```

Important outputs:

- `fold0_train_cases.txt`
- `fold0_validation_cases.txt`
- `kalasanty_valid_cases.txt`
- `summary.json`
- `summary.csv`

Kalasanty fold IDs are base PDB IDs. This script maps those base IDs back to our one-H5-per-scPDB-entry cache files, while keeping fold separation by base PDB ID.
