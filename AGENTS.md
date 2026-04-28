# AGENTS.md

## Purpose
This repository generates HDF5 cache files for protein-ligand binding-site experiments. Treat it as the preprocessing/data-engineering repo for PDBBind and scPDB inputs.

## Active Pipeline
- Prefer `generate_cache_pdbbind_resolution_fix.py` for current cache generation unless the user explicitly asks for an older script.
- Keep the cache schema compatible with `3dunet-apbs/3dunet_configurable/dataset.py`: features under `features/` and labels under `label/`.
- Generated HDF5 files, APBS temporary directories, `.dx`, `.pqr`, local conda environments, and visualization outputs are artifacts. Do not commit them.

## Data Integrity Rules
- Do not use ligand-derived information as a model input for de novo binding-site prediction unless the experiment is explicitly ligand-conditioned.
- If adding a new feature, document whether it is available at inference time.
- Keep physical units explicit. For resolution-fix caches, the physical extent is fixed at 161 Angstrom and voxel size is `161.0 / box_size`.
- Avoid silently changing label definitions. If label generation changes, write a new label name or add a clear config/script note.

## Implementation Rules
- Use structured parsers/APIs for PDB, MOL2, HDF5, and YAML instead of ad hoc string parsing where practical.
- Keep dataset paths configurable when adding new scripts. Do not hard-code new machine-specific paths unless matching the existing local experiment style and called out clearly.
- Preserve failed-case logging. Cache generation should skip existing successful outputs and remove partial `.h5` files on failure.
- Do not run large cache generation jobs unless the user explicitly asks.

## Verification
- For code-only changes, run at least an import/syntax check for touched scripts.
- For cache schema changes, inspect one generated `.h5` file and report feature names, label names, tensor shapes, positive voxel counts, and min/max statistics for continuous features.
