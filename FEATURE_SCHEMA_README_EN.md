# Deep-APBS H5 Feature Schema

This document describes the HDF5 feature schema used by the Deep-APBS cache files. The goal is to make the dataset self-describing, reusable, and safe for protein-only blind binding-site prediction experiments.

## H5 File Layout

| Group | Meaning | Should it be used as model input? |
|---|---|---|
| `features/` | Protein-derived candidate input channels for blind binding-site prediction. | Yes, depending on the selected ablation. |
| `label/` | Supervision and evaluation targets. | No, target only. |
| `auxiliary/` | Ligand-derived or debug channels kept for visualization and sanity checks. | No. These channels leak ligand information. |
| root attributes | Grid, APBS, and schema metadata. | No. |

`features/electrostatic_grid` is a legacy name and is not part of the public schema. It is used only as an input to create the explicit `electrostatic_grid_v1_ligand_proximal_chains_7A_raw` channel, then removed.

## Grid and Physical Scale

| Setting | Meaning | Approximate resolution |
|---|---|---|
| `box36_span70` | Small Kalasanty/PUResNet-like box with 70 A physical span. | `70 / 35 = 2.00 A/voxel` |
| `box72_span120` | Wider context while preserving better resolution than box36. | `120 / 71 = 1.69 A/voxel` |
| `box161_span160` | Near-native APBS field resolution. | `160 / 160 = 1.00 A/voxel` |

APBS is computed on a 161-point grid spanning 160 A. Smaller model grids resample this APBS field onto the target grid.

## APBS Source Definitions

| Prefix | Definition | Role |
|---|---|---|
| `electrostatic_grid_v1_ligand_proximal_chains_7A_*` | Legacy APBS computed after selecting ligand-proximal protein chains/residues. | Reproduces and compares against legacy experiments. |
| `electrostatic_grid_v2_full_protein_*` | APBS computed from the full protein structure without ligand-proximal chain trimming. | Cleaner representation for protein-only blind prediction. |

## APBS Representations

The same APBS potential field is stored in multiple normalized representations because electrostatic values may contain strong outliers and neural networks can be sensitive to input scale.

| Suffix | Description | Value range | Notes |
|---|---|---|---|
| `raw` | Raw APBS electrostatic potential resampled onto the target grid. | Unbounded | Kept for analysis and alternative normalization. |
| `clip5_minmax` | Values clipped to `[-5, 5]`, then scaled to `[0, 1]`. | `[0, 1]` | Emphasizes weak to moderate potential regions. |
| `clip10_minmax` | Values clipped to `[-10, 10]`, then scaled to `[0, 1]`. | `[0, 1]` | Preserves moderate electrostatic signal. |
| `clip20_minmax` | Values clipped to `[-20, 20]`, then scaled to `[0, 1]`. | `[0, 1]` | Preserves stronger potential differences. |
| `clip20_signed` | Values clipped to `[-20, 20]`, then scaled to `[-1, 1]`. | `[-1, 1]` | Preserves positive and negative sign. |
| `full_signed150` | Raw values divided by `150` without clipping. | Usually near `[-1, 1]`, theoretically unbounded | Keeps high-magnitude information instead of removing it. |
| `clip150_signed` | Values clipped to `[-150, 150]`, then scaled to `[-1, 1]`. | `[-1, 1]` | Bounded version of the full signed representation. |
| `positive_clip20` | Positive component clipped to `[0, 20]`, then scaled to `[0, 1]`. | `[0, 1]` | Separates positive electrostatic regions. |
| `negative_clip20` | Negative component stored as positive magnitude, clipped to `[0, 20]`, then scaled to `[0, 1]`. | `[0, 1]` | Separates negative electrostatic regions. |
| `gradient_magnitude_robust` | Spatial APBS gradient magnitude normalized to `[0, 1]` with a robust percentile scale. | `[0, 1]` | Highlights where electrostatic potential changes rapidly. |
| `clip20_signed_surface_weighted` | `clip20_signed * protein_proximity_exp3`. | Around `[-1, 1]` | Emphasizes near-protein electrostatic signal. |
| `full_signed150_surface_weighted` | `full_signed150 * protein_proximity_exp3`. | Usually around `[-1, 1]` | Emphasizes near-protein signed APBS signal without clipping. |

## Public Model Features

### Atomic Channels

| Feature | Description |
|---|---|
| `atomic_B` | Voxelized boron atom occupancy. |
| `atomic_C` | Voxelized carbon atom occupancy. |
| `atomic_N` | Voxelized nitrogen atom occupancy. |
| `atomic_O` | Voxelized oxygen atom occupancy. |
| `atomic_P` | Voxelized phosphorus atom occupancy. |
| `atomic_S` | Voxelized sulfur atom occupancy. |
| `atomic_Se` | Voxelized selenium atom occupancy. |
| `atomic_acceptor` | Hydrogen-bond acceptor atom/property channel. |
| `atomic_aromatic` | Aromatic atom or aromatic ring participation channel. |
| `atomic_donor` | Hydrogen-bond donor atom/property channel. |
| `atomic_halogen` | Halogen atom/property channel. |
| `atomic_heavydegree` | Heavy-atom bonding degree channel. |
| `atomic_heterodegree` | Heteroatom bonding degree channel. |
| `atomic_hyb` | Atom hybridization descriptor channel. |
| `atomic_hydrophobic` | Atom-level hydrophobic property channel. |
| `atomic_metal` | Metal atom/property channel. |
| `atomic_molcode` | Molecule-code channel emitted by the atom featurizer. |
| `atomic_partialcharge` | Atom partial charge channel. |
| `atomic_ring` | Ring-membership descriptor channel. |

### Geometry and Physicochemical Channels

| Feature | Description |
|---|---|
| `shape` | Binary protein occupancy mask. It represents the protein volume and overall shape. |
| `hydrophobicity` | Residue hydrophobicity values rasterized from protein atoms onto the grid. |
| `dist_to_surface` | Distance from each grid point to the nearest protein atom or surface proxy. It can encode pocket depth and surface proximity. |
| `protein_proximity_exp3` | Protein-proximity channel computed as `exp(-dist_to_surface / 3A)`. It softly emphasizes near-protein regions. |
| `protein_near_shell_0_3A` | Binary near-protein shell where `dist_to_surface` is between 0 and 3 A. |
| `protein_near_shell_3_6A` | Binary outer near-protein shell where `dist_to_surface` is greater than 3 A and up to 6 A. |
| `hydrophobicity_surface_weighted` | `hydrophobicity * protein_proximity_exp3`. It emphasizes hydrophobic signal near the protein/surface. |

### APBS v1 Channels

All channels below are derived from the ligand-proximal legacy APBS source:

- `electrostatic_grid_v1_ligand_proximal_chains_7A_raw`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_clip5_minmax`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_clip10_minmax`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_clip20_minmax`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_clip20_signed`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_full_signed150`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_clip150_signed`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_positive_clip20`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_negative_clip20`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_gradient_magnitude_robust`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_clip20_signed_surface_weighted`
- `electrostatic_grid_v1_ligand_proximal_chains_7A_full_signed150_surface_weighted`

### APBS v2 Channels

All channels below are derived from full-protein APBS:

- `electrostatic_grid_v2_full_protein_raw`
- `electrostatic_grid_v2_full_protein_clip5_minmax`
- `electrostatic_grid_v2_full_protein_clip10_minmax`
- `electrostatic_grid_v2_full_protein_clip20_minmax`
- `electrostatic_grid_v2_full_protein_clip20_signed`
- `electrostatic_grid_v2_full_protein_full_signed150`
- `electrostatic_grid_v2_full_protein_clip150_signed`
- `electrostatic_grid_v2_full_protein_positive_clip20`
- `electrostatic_grid_v2_full_protein_negative_clip20`
- `electrostatic_grid_v2_full_protein_gradient_magnitude_robust`
- `electrostatic_grid_v2_full_protein_clip20_signed_surface_weighted`
- `electrostatic_grid_v2_full_protein_full_signed150_surface_weighted`

## Labels

| Label | Description | Use |
|---|---|---|
| `binding_site_calculated` | Binding-site label calculated from protein-ligand proximity by the cache pipeline. | Supervision or evaluation. |
| `binding_site_in_dataset` | Binding-site label provided by the source dataset when available. | Supervision or evaluation. |
| `binding_site_fpocket_selected` | Selected fpocket pocket label for external benchmark datasets. | External benchmark evaluation. |

## Auxiliary Channels

| Auxiliary channel | Description | Use as model input? |
|---|---|---|
| `ligand` | Voxelized ligand mask. | No. Visualization/debug only. |
| `dist_to_ligand` | Distance from each grid point to ligand atoms. | No. It directly leaks target information. |

## Recommended Ablation Groups

| Group | Feature set |
|---|---|
| Shape only | `shape` |
| APBS only | One APBS representation, preferably the strongest candidate from current experiments such as `electrostatic_grid_v2_full_protein_full_signed150` or the best validated v1/v2 representation. |
| Shape + APBS | `shape` + selected APBS channel. |
| Selected chemistry | `atomic_C`, `atomic_N`, `atomic_O`, `atomic_S`, `atomic_acceptor`, `atomic_donor`, `atomic_aromatic`, `atomic_hydrophobic`, `atomic_partialcharge`, `hydrophobicity` |
| APBS + selected chemistry | Selected APBS channel + selected chemistry. |
| APBS + shape + selected chemistry | Selected APBS channel + `shape` + selected chemistry. |
| Surface/hydro add-on | Strong groups above plus `dist_to_surface` and `hydrophobicity`. |

## Leakage Rule

In the intended blind prediction setting, the user provides only the protein structure. Therefore, `label/`, `auxiliary/ligand`, and `auxiliary/dist_to_ligand` must never be used as input features. They are kept only for validation, visualization, and error analysis.
