# AGENTS.md

## Purpose
This repository generates HDF5 cache files for protein-ligand binding-site experiments. Treat it as the preprocessing/data-engineering repo for PDBBind and scPDB inputs.

## Communication
- Answer the user in the language they use. Keep code, identifiers, comments, config keys, logs, and commit-style technical text in English.

## Active Pipeline
- Prefer `generate_cache_pdbbind_gridfix.py` for new PDBBind box-size scaling experiments, especially local `box72` APBS-downsampled caches.
- Keep `generate_cache_pdbbind_resolution_fix.py` only as the legacy resolution-fix path for reproducing older cache outputs.
- Before broad gridfix generation, use `audit_gridfix_coverage.py` to produce/review the case list for the intended physical span.
- Keep the cache schema compatible with `3dunet-apbs/3dunet_configurable/dataset.py`: features under `features/` and labels under `label/`.
- Generated HDF5 files, APBS temporary directories, `.dx`, `.pqr`, local conda environments, and visualization outputs are artifacts. Do not commit them.

## Project Direction
- Current priority is to finish the planned experimental work packages and produce reliable caches/results for the thesis before doing framework/product work.
- After the experiment results are stable, support turning the broader project into a protein binding-site prediction framework inspired by `pytorch-3dunet`: config-driven cache generation, training, prediction, evaluation, standardized HDF5/cache format, pretrained checkpoints, and visualization.
- The thesis does not have to depend on the framework being complete. Treat framework readiness as a follow-up deliverable and CV-strengthening artifact unless the user explicitly changes priorities.
- Even before the framework phase, write cache scripts in a framework-compatible style: clear CLI arguments, reproducible configs, stable output folders, machine-readable logs/summaries, and no hidden local assumptions.
- Preserve legacy features while adding improved APBS/MOL2/atomic feature paths so old experiments remain reproducible and new ablations are clean.

## Data Integrity Rules
- Do not use ligand-derived information as a model input for de novo binding-site prediction unless the experiment is explicitly ligand-conditioned.
- If adding a new feature, document whether it is available at inference time.
- Keep physical units explicit. For gridfix caches, APBS is computed as 161 grid points spanning 160 Angstrom, and target voxel size is `160.0 / (box_size - 1)`.
- Avoid silently changing label definitions. If label generation changes, write a new label name or add a clear config/script note.

## Implementation Rules
- Use structured parsers/APIs for PDB, MOL2, HDF5, and YAML instead of ad hoc string parsing where practical.
- Keep dataset paths configurable when adding new scripts. Do not hard-code new machine-specific paths unless matching the existing local experiment style and called out clearly.
- Preserve failed-case logging. Cache generation should skip existing successful outputs and remove partial `.h5` files on failure.
- Do not run large cache generation jobs unless the user explicitly asks.

## Verification
- For code-only changes, run at least an import/syntax check for touched scripts.
- For cache schema changes, inspect one generated `.h5` file and report feature names, label names, tensor shapes, positive voxel counts, and min/max statistics for continuous features.
