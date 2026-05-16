"""Feature schema metadata for Deep-APBS H5 cache files."""

from __future__ import annotations

FEATURE_SCHEMA_VERSION = "deep_apbs_feature_schema_v2"

FILE_ATTRS = {
    "feature_schema_version": FEATURE_SCHEMA_VERSION,
    "feature_schema_summary": (
        "features/ contains model-input candidate channels; auxiliary/ contains "
        "debug or leakage-prone channels that must not be used as blind prediction inputs; "
        "label/ contains supervision and evaluation targets."
    ),
    "features_are_model_input_candidates": "true",
    "auxiliary_is_not_for_training": "true",
    "labels_are_ground_truth_targets": "true",
    "apbs_v1_definition": (
        "APBS electrostatic potential computed on the ligand-proximal selected protein "
        "chains/residues using the legacy 7A ligand-proximity selection pipeline."
    ),
    "apbs_v2_definition": "APBS electrostatic potential computed on the full protein structure.",
    "legacy_electrostatic_grid_policy": (
        "features/electrostatic_grid is intentionally removed after explicit v1 APBS "
        "channels are created."
    ),
}

ATOMIC_DESCRIPTIONS = {
    "atomic_B": "Voxelized boron atom occupancy channel.",
    "atomic_C": "Voxelized carbon atom occupancy channel.",
    "atomic_N": "Voxelized nitrogen atom occupancy channel.",
    "atomic_O": "Voxelized oxygen atom occupancy channel.",
    "atomic_P": "Voxelized phosphorus atom occupancy channel.",
    "atomic_S": "Voxelized sulfur atom occupancy channel.",
    "atomic_Se": "Voxelized selenium atom occupancy channel.",
    "atomic_acceptor": "Voxelized hydrogen-bond acceptor property from atom-level featurization.",
    "atomic_aromatic": "Voxelized aromatic atom/ring participation property.",
    "atomic_donor": "Voxelized hydrogen-bond donor property from atom-level featurization.",
    "atomic_halogen": "Voxelized halogen atom/property channel.",
    "atomic_heavydegree": "Voxelized heavy-atom degree, describing local heavy-atom bonding count.",
    "atomic_heterodegree": "Voxelized heteroatom degree, describing local heteroatom bonding count.",
    "atomic_hyb": "Voxelized atom hybridization descriptor.",
    "atomic_hydrophobic": "Voxelized atom-level hydrophobic property.",
    "atomic_metal": "Voxelized metal atom/property channel.",
    "atomic_molcode": "Voxelized molecule-code channel emitted by the atom featurizer.",
    "atomic_partialcharge": "Voxelized atom partial charge descriptor.",
    "atomic_ring": "Voxelized ring-membership descriptor.",
}

GENERAL_FEATURES = {
    "shape": {
        "category": "geometry",
        "description": "Binary protein occupancy mask on the grid.",
        "recommended_use": "model_input",
    },
    "hydrophobicity": {
        "category": "physicochemical",
        "description": "Residue hydrophobicity values rasterized from protein atoms onto the grid.",
        "recommended_use": "model_input",
    },
    "dist_to_surface": {
        "category": "geometry",
        "description": "Distance from each grid point to the nearest protein atom/surface proxy.",
        "recommended_use": "model_input_candidate",
    },
    "protein_proximity_exp3": {
        "category": "geometry_derived",
        "description": "Exponential protein-proximity channel computed as exp(-dist_to_surface / 3A).",
        "recommended_use": "model_input_candidate",
    },
    "protein_near_shell_0_3A": {
        "category": "geometry_derived",
        "description": "Binary near-protein shell where dist_to_surface is between 0A and 3A.",
        "recommended_use": "model_input_candidate",
    },
    "protein_near_shell_3_6A": {
        "category": "geometry_derived",
        "description": "Binary near-protein shell where dist_to_surface is greater than 3A and up to 6A.",
        "recommended_use": "model_input_candidate",
    },
    "hydrophobicity_surface_weighted": {
        "category": "physicochemical_derived",
        "description": "Hydrophobicity weighted by protein_proximity_exp3 to emphasize near-surface hydrophobic signal.",
        "recommended_use": "model_input_candidate",
    },
}

APBS_REPRESENTATIONS = {
    "raw": "Raw APBS electrostatic potential resampled onto the target grid.",
    "clip5_minmax": "APBS values clipped to [-5, 5] and min-max scaled to [0, 1].",
    "clip10_minmax": "APBS values clipped to [-10, 10] and min-max scaled to [0, 1].",
    "clip20_minmax": "APBS values clipped to [-20, 20] and min-max scaled to [0, 1].",
    "clip20_signed": "APBS values clipped to [-20, 20] and scaled to [-1, 1].",
    "full_signed150": "Unclipped APBS values divided by 150, preserving signed magnitude.",
    "clip150_signed": "APBS values clipped to [-150, 150] and scaled to [-1, 1].",
    "positive_clip20": "Positive APBS component clipped to [0, 20] and scaled to [0, 1].",
    "negative_clip20": "Negative APBS component stored as positive magnitude, clipped to [0, 20] and scaled to [0, 1].",
    "gradient_magnitude_robust": "Robustly normalized magnitude of the spatial APBS gradient, encoding how quickly electrostatic potential changes over the grid.",
    "clip20_signed_surface_weighted": "clip20_signed APBS multiplied by protein_proximity_exp3 to emphasize near-protein electrostatic signal.",
    "full_signed150_surface_weighted": "full_signed150 APBS multiplied by protein_proximity_exp3 to emphasize near-protein signed electrostatic signal.",
}

LABEL_DESCRIPTIONS = {
    "binding_site_calculated": {
        "category": "label",
        "description": "Binding-site label calculated from protein-ligand proximity by the cache pipeline.",
        "recommended_use": "supervision_or_evaluation",
    },
    "binding_site_in_dataset": {
        "category": "label",
        "description": "Binding-site label provided by the source dataset when available.",
        "recommended_use": "supervision_or_evaluation",
    },
    "binding_site_fpocket_selected": {
        "category": "label",
        "description": "Selected fpocket pocket label for external benchmark datasets.",
        "recommended_use": "external_benchmark_evaluation",
    },
}

AUXILIARY_DESCRIPTIONS = {
    "ligand": {
        "category": "auxiliary_leakage",
        "description": "Voxelized ligand mask. Kept only for visualization/debugging; not available during blind protein-only prediction.",
        "recommended_use": "do_not_use_for_training",
    },
    "dist_to_ligand": {
        "category": "auxiliary_leakage",
        "description": "Distance from each grid point to ligand atoms. Kept only for visualization/debugging; it leaks target information.",
        "recommended_use": "do_not_use_for_training",
    },
    "electrostatic_grid": {
        "category": "legacy",
        "description": "Deprecated legacy raw APBS channel. Public schema removes this in favor of explicit v1 APBS names.",
        "recommended_use": "deprecated",
    },
}


def apbs_feature_metadata(name: str) -> dict[str, str] | None:
    prefixes = {
        "electrostatic_grid_v1_ligand_proximal_chains_7A_": (
            "apbs_v1_ligand_proximal",
            FILE_ATTRS["apbs_v1_definition"],
        ),
        "electrostatic_grid_v2_full_protein_": (
            "apbs_v2_full_protein",
            FILE_ATTRS["apbs_v2_definition"],
        ),
    }
    for prefix, (category, source_definition) in prefixes.items():
        if name.startswith(prefix):
            representation = name[len(prefix) :]
            description = APBS_REPRESENTATIONS.get(representation)
            if description is None:
                return None
            return {
                "category": category,
                "description": f"{description} Source: {source_definition}",
                "recommended_use": "model_input_candidate",
                "apbs_source_definition": source_definition,
                "apbs_representation": representation,
            }
    return None


def feature_metadata(name: str) -> dict[str, str] | None:
    if name in ATOMIC_DESCRIPTIONS:
        return {
            "category": "atomic",
            "description": ATOMIC_DESCRIPTIONS[name],
            "recommended_use": "model_input_candidate",
        }
    if name in GENERAL_FEATURES:
        return GENERAL_FEATURES[name]
    if name in AUXILIARY_DESCRIPTIONS:
        return AUXILIARY_DESCRIPTIONS[name]
    return apbs_feature_metadata(name)


def apply_dataset_metadata(dataset, metadata: dict[str, str]) -> None:
    for key, value in metadata.items():
        dataset.attrs[key] = str(value)


def apply_feature_schema_metadata(h5f) -> None:
    for key, value in FILE_ATTRS.items():
        h5f.attrs[key] = str(value)

    if "features" in h5f:
        for name, dataset in h5f["features"].items():
            metadata = feature_metadata(name)
            if metadata is not None:
                apply_dataset_metadata(dataset, metadata)

    if "auxiliary" in h5f:
        for name, dataset in h5f["auxiliary"].items():
            metadata = AUXILIARY_DESCRIPTIONS.get(name)
            if metadata is not None:
                apply_dataset_metadata(dataset, metadata)

    if "label" in h5f:
        for name, dataset in h5f["label"].items():
            metadata = LABEL_DESCRIPTIONS.get(name)
            if metadata is not None:
                apply_dataset_metadata(dataset, metadata)
