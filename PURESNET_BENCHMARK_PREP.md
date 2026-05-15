# PUResNet Coach420 and BU48 Benchmark Preparation

This note documents how the local Coach420/BU48 benchmark staging was prepared.

## Sources

- PUResNet repository benchmark zips:
  - `coach.zip`
  - `BU48.zip`
- P2Rank dataset repository:
  - `rdk/p2rank-datasets`
  - Provides `coach420-fpocket.ds`, `joined-fpocket.ds`, and fpocket output folders.

The P2Rank dataset README states that stored fpocket predictions were produced
with fpocket v1.0 using default parameters.

## Why P2Rank fpocket Outputs Are Used

PUResNet's public benchmark zips contain protein-ligand pairs, but they do not
contain explicit ground-truth pocket volumes. In correspondence, the PUResNet
author described the benchmark evaluation approach as:

- DCC can be calculated between the predicted pocket center and ligand center.
- DVO can use fpocket pockets.
- The relevant fpocket pocket can be selected by DCA, i.e. choosing the pocket
  associated with the ligand.
- The author pointed to `rdk/p2rank-datasets`, where fpocket results are already
  provided.

The PUResNet paper reports Coach420/BU48 results, but the exact fpocket selection
procedure is not fully described in the manuscript text. Therefore this local
preparation keeps the source and selection method explicit in `manifest.csv`.

## Prepared Subsets

The preparation script keeps the PUResNet public benchmark subsets:

- `coach420_puresnet`: 298 cases from PUResNet `coach.zip`
- `bu48_puresnet`: 62 cases from PUResNet `BU48.zip`

The P2Rank repository contains more Coach420 and BU48 entries than these subsets.
Those extra entries are not staged here, because the goal is to reproduce the
PUResNet benchmark subset first.

## Command

Run from the `generate_cache` repository:

```bash
.conda/bin/python prepare_puresnet_benchmarks.py \
  --output-root /Users/tevfik/Sandbox/github/PHD/data/external_benchmarks/puresnet_prepared \
  --overwrite
```

## Output

Root:

```text
/Users/tevfik/Sandbox/github/PHD/data/external_benchmarks/puresnet_prepared
```

Files:

- `manifest.csv`: one row per prepared benchmark case.
- `coach420_puresnet_cases.txt`: 298 prepared Coach case ids.
- `bu48_puresnet_cases.txt`: 62 prepared BU48 case ids.
- `README.md`: short artifact-level description.

Per-case files:

```text
<output-root>/<dataset>/<case_id>/protein.pdb
<output-root>/<dataset>/<case_id>/ligand.pdb
<output-root>/<dataset>/<case_id>/fpocket_selected_pocket_atoms.pdb
<output-root>/<dataset>/<case_id>/fpocket_selected_pocket_vertices.pqr
```

BU48 P2Rank fpocket folders generally do not include `pocket*_vert.pqr` files,
so `fpocket_selected_pocket_vertices.pqr` is absent for those cases.

## Pocket Selection

Default selection method:

```text
center_dca
```

For every fpocket pocket:

1. Read `pocket*_atm.pdb`.
2. Compute the centroid of the pocket atoms.
3. Compute the distance from that centroid to the ligand centroid.
4. Select the pocket with the smallest distance.

The manifest also records the minimum atom-to-atom ligand distance as a
diagnostic, but this is not the default selection score.

## Verification

Last verified counts:

```text
coach420_puresnet: 298 ok, 0 failed
bu48_puresnet: 62 ok, 0 failed
total: 360 ok, 0 failed
```

The prepared files are benchmark staging inputs. They are not HDF5 training
caches yet. A later cache-generation step should convert the selected fpocket
pocket into an explicit label/mask if DVO or voxel-level evaluation requires it.
