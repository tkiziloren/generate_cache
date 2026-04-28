import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
from plotly.subplots import make_subplots
import plotly.graph_objects as go



# Set which features are main (isosurface/grid) and which are atomic (point)
MAIN_FEATURES = ["electrostatic_grid", "dist_to_surface", "dist_to_ligand", "shape"]
ATOM_FEATURES = [
    "atomic_B", "atomic_C", "atomic_N", "atomic_O", "atomic_P", "atomic_S", "atomic_Se",
    "atomic_acceptor", "atomic_aromatic", "atomic_donor", "atomic_halogen", "atomic_heavydegree",
    "atomic_heterodegree", "atomic_hyb", "atomic_hydrophobic", "atomic_metal", "atomic_molcode",
    "atomic_partialcharge", "atomic_ring"
]
import math
def plot_all_feature_histograms(h5_path, mask_type='surface'):
    with h5py.File(h5_path, "r") as h5f:
        # Feature listesi
        features = list(h5f['features'].keys())
        label_names = ["binding_site_calculated", "binding_site_in_dataset"]
        total_plots = len(features) + len(label_names)

        # Grid boyutları
        cols = 6
        rows = math.ceil((total_plots + 2) / cols)
        plt.figure(figsize=(4*cols, 3*rows))

        # Mask oluştur
        shape_mask = h5f["features/shape"][:]
        dist2surface = h5f["features/dist_to_surface"][:]
        if mask_type == 'surface':
            mask = (dist2surface < 2.0)
        elif mask_type == 'protein':
            mask = (shape_mask > 0)
        else:
            mask = np.ones_like(shape_mask, dtype=bool)

        # Önce feature histogramları
        for i, feat in enumerate(features):
            ax = plt.subplot(rows, cols, i+1)
            data = h5f['features'][feat][:]
            ax.hist(data[mask].flatten(), bins=40, alpha=0.7)
            ax.set_title(feat)

        # Ardından label histogramları
        for j, lbl in enumerate(label_names, start=len(features)):
            ax = plt.subplot(rows, cols, j+1)
            data_lbl = h5f[f"label/{lbl}"][:].flatten()
            ax.hist(data_lbl, bins=40, alpha=0.7)
            ax.set_title(lbl)

        plt.tight_layout()
        plt.show()



def plot_all_features_in_grid(
    h5_path, ncols=6, atomic_marker_size=4,
    apbs_min=-15, apbs_max=15,  # APBS aralığını elle belirle!
    percentile=98,
    atomic_prefix="atomic_"
):
    # --- Key features ve atomic features ayrımı ---
    key_features = [
        "apbs", "dist_to_surface", "dist_to_ligand", "shape", "protein surface binding with calculated binding site", "protein surface binding with provided binding site"
    ]
    atomic_features = []
    with h5py.File(h5_path, "r") as h5f:
        for feat in h5f['features']:
            if feat.startswith(atomic_prefix):  # sadece atomic_ ile başlayanları ekle
                atomic_features.append(feat)
    n_key = len(key_features)
    n_atomic = len(atomic_features)
    total = n_key + n_atomic
    nrows = 1 + int(np.ceil(n_atomic / ncols))

    fig = make_subplots(
        rows=nrows, cols=ncols,
        specs=[[{'type': 'scene'}]*ncols for _ in range(nrows)],
        subplot_titles=key_features + atomic_features,
        horizontal_spacing=0.02,
        vertical_spacing=0.03
    )

    with h5py.File(h5_path, 'r') as h5f:
        # --- 1. APBS
        idx = 0
        row, col = 1, 1
        apbs = h5f['features/electrostatic_grid'][:]
        fig.add_trace(go.Isosurface(
            x=np.arange(apbs.shape[0]).repeat(apbs.shape[1]*apbs.shape[2]),
            y=np.tile(np.arange(apbs.shape[1]).repeat(apbs.shape[2]), apbs.shape[0]),
            z=np.tile(np.arange(apbs.shape[2]), apbs.shape[0]*apbs.shape[1]),
            value=apbs.flatten(),
            isomin=apbs_min, isomax=apbs_max,
            surface_count=2,
            colorscale="RdBu",
            opacity=0.55,
            showscale=True,
            colorbar=dict(
                title=dict(text="APBS", font=dict(size=18)),
                thickness=14,
                len=0.7,
                x=1.06,
                xanchor="left",
                outlinewidth=1,
                tickvals=[apbs_min, 0, apbs_max],
                ticktext=[f"{apbs_min}", "0", f"{apbs_max}"]
            ),
            caps=dict(x_show=False, y_show=False, z_show=False),
            name="apbs",
            showlegend=False
        ), row=row, col=col)

        # --- 2. dist2surface
        idx += 1
        row, col = 1, 2
        dist2surface = h5f['features/dist_to_surface'][:]
        dist2surface_nonzero = dist2surface[dist2surface != 0]
        d2s_isovalue = np.percentile(dist2surface_nonzero, percentile) if dist2surface_nonzero.size > 0 else 0.1
        fig.add_trace(go.Isosurface(
            x=np.arange(dist2surface.shape[0]).repeat(dist2surface.shape[1]*dist2surface.shape[2]),
            y=np.tile(np.arange(dist2surface.shape[1]).repeat(dist2surface.shape[2]), dist2surface.shape[0]),
            z=np.tile(np.arange(dist2surface.shape[2]), dist2surface.shape[0]*dist2surface.shape[1]),
            value=dist2surface.flatten(),
            isomin=d2s_isovalue, isomax=dist2surface.max(),
            surface_count=2,
            colorscale="Viridis",
            opacity=0.7,
            showscale=False,
            caps=dict(x_show=False, y_show=False, z_show=False),
            name="dist2surface",
            showlegend=False
        ), row=row, col=col)

        # --- 3. dist2ligand
        idx += 1
        row, col = 1, 3
        dist2ligand = h5f['features/dist_to_ligand'][:]
        dist2ligand_nonzero = dist2ligand[dist2ligand != 0]
        d2l_isovalue = np.percentile(dist2ligand_nonzero, percentile) if dist2ligand_nonzero.size > 0 else 0.1
        fig.add_trace(go.Isosurface(
            x=np.arange(dist2ligand.shape[0]).repeat(dist2ligand.shape[1]*dist2ligand.shape[2]),
            y=np.tile(np.arange(dist2ligand.shape[1]).repeat(dist2ligand.shape[2]), dist2ligand.shape[0]),
            z=np.tile(np.arange(dist2ligand.shape[2]), dist2ligand.shape[0]*dist2ligand.shape[1]),
            value=dist2ligand.flatten(),
            isomin=d2l_isovalue, isomax=dist2ligand.max(),
            surface_count=2,
            colorscale="Magma",
            opacity=0.7,
            showscale=False,
            caps=dict(x_show=False, y_show=False, z_show=False),
            name="dist2ligand",
            showlegend=False
        ), row=row, col=col)

        # --- 4. shape
        idx += 1
        row, col = 1, 4
        shape = h5f['features/shape'][:]
        fig.add_trace(go.Isosurface(
            x=np.arange(shape.shape[0]).repeat(shape.shape[1]*shape.shape[2]),
            y=np.tile(np.arange(shape.shape[1]).repeat(shape.shape[2]), shape.shape[0]),
            z=np.tile(np.arange(shape.shape[2]), shape.shape[0]*shape.shape[1]),
            value=shape.flatten(),
            isomin=0.5, isomax=1.0,
            opacity=0.18,
            surface_count=1,
            colorscale="Blues",
            showscale=False,
            caps=dict(x_show=False, y_show=False, z_show=False),
            name="shape",
            showlegend=False
        ), row=row, col=col)

        # --- 5. protein_surface_binding: shape, binding site ve ligand aynı panelde
        idx += 1
        row, col = 1, 5
        label = h5f["label/binding_site_calculated"][:]
        label_provided = h5f["label/binding_site_in_dataset"][:]
        ligand_mask = h5f["features/ligand"][:]

        # Shape (Protein surface)
        fig.add_trace(go.Isosurface(
            x=np.arange(shape.shape[0]).repeat(shape.shape[1]*shape.shape[2]),
            y=np.tile(np.arange(shape.shape[1]).repeat(shape.shape[2]), shape.shape[0]),
            z=np.tile(np.arange(shape.shape[2]), shape.shape[0]*shape.shape[1]),
            value=shape.flatten(),
            isomin=0.5, isomax=1.0,
            opacity=0.15,
            surface_count=1,
            colorscale="Blues",
            showscale=False,
            caps=dict(x_show=False, y_show=False, z_show=False),
            name='protein_surface',
            showlegend=False
        ), row=row, col=col)
        # Binding site
        fig.add_trace(go.Isosurface(
            x=np.arange(label.shape[0]).repeat(label.shape[1]*label.shape[2]),
            y=np.tile(np.arange(label.shape[1]).repeat(label.shape[2]), label.shape[0]),
            z=np.tile(np.arange(label.shape[2]), label.shape[0]*label.shape[1]),
            value=label.flatten(),
            isomin=0.5, isomax=1.0,
            opacity=0.54,
            surface_count=1,
            colorscale="Reds",
            showscale=False,
            caps=dict(x_show=False, y_show=False, z_show=False),
            name='binding_site',
            showlegend=False
        ), row=row, col=col)
        # Ligand (aynı panelde, renk: yeşil)
        fig.add_trace(go.Isosurface(
            x=np.arange(ligand_mask.shape[0]).repeat(ligand_mask.shape[1]*ligand_mask.shape[2]),
            y=np.tile(np.arange(ligand_mask.shape[1]).repeat(ligand_mask.shape[2]), ligand_mask.shape[0]),
            z=np.tile(np.arange(ligand_mask.shape[2]), ligand_mask.shape[0]*ligand_mask.shape[1]),
            value=ligand_mask.flatten(),
            isomin=0.5, isomax=1.0,
            opacity=0.68,
            surface_count=1,
            colorscale="Greens",
            showscale=False,
            caps=dict(x_show=False, y_show=False, z_show=False),
            name='ligand',
            showlegend=False
        ), row=row, col=col)

        # --- 6. Ligand (Ayrı panel) ---
        idx += 1
        row, col = 1, 6
        
        # Shape (Protein surface)
        fig.add_trace(go.Isosurface(
            x=np.arange(shape.shape[0]).repeat(shape.shape[1]*shape.shape[2]),
            y=np.tile(np.arange(shape.shape[1]).repeat(shape.shape[2]), shape.shape[0]),
            z=np.tile(np.arange(shape.shape[2]), shape.shape[0]*shape.shape[1]),
            value=shape.flatten(),
            isomin=0.5, isomax=1.0,
            opacity=0.15,
            surface_count=1,
            colorscale="Blues",
            showscale=False,
            caps=dict(x_show=False, y_show=False, z_show=False),
            name='protein_surface',
            showlegend=False
        ), row=row, col=col)
        # Binding site
        fig.add_trace(go.Isosurface(
            x=np.arange(label_provided.shape[0]).repeat(label_provided.shape[1]*label_provided.shape[2]),
            y=np.tile(np.arange(label_provided.shape[1]).repeat(label_provided.shape[2]), label_provided.shape[0]),
            z=np.tile(np.arange(label_provided.shape[2]), label_provided.shape[0]*label_provided.shape[1]),
            value=label_provided.flatten(),
            isomin=0.5, isomax=1.0,
            opacity=0.54,
            surface_count=1,
            colorscale="Reds",
            showscale=False,
            caps=dict(x_show=False, y_show=False, z_show=False),
            name='binding_site',
            showlegend=False
        ), row=row, col=col)
        # Ligand (aynı panelde, renk: yeşil)
        fig.add_trace(go.Isosurface(
            x=np.arange(ligand_mask.shape[0]).repeat(ligand_mask.shape[1]*ligand_mask.shape[2]),
            y=np.tile(np.arange(ligand_mask.shape[1]).repeat(ligand_mask.shape[2]), ligand_mask.shape[0]),
            z=np.tile(np.arange(ligand_mask.shape[2]), ligand_mask.shape[0]*ligand_mask.shape[1]),
            value=ligand_mask.flatten(),
            isomin=0.5, isomax=1.0,
            opacity=0.68,
            surface_count=1,
            colorscale="Greens",
            showscale=False,
            caps=dict(x_show=False, y_show=False, z_show=False),
            name='ligand',
            showlegend=False
        ), row=row, col=col)

        # --- Atomic Features ---
        for j, feat in enumerate(atomic_features):
            idx = n_key + j
            row = idx // ncols + 1
            col = idx % ncols + 1
            data = h5f[f'features/{feat}'][:]
            pts = np.where(data > 0.5)
            if pts[0].size > 0:
                fig.add_trace(go.Scatter3d(
                    x=pts[0], y=pts[1], z=pts[2],
                    mode="markers",
                    marker=dict(size=atomic_marker_size, color="purple", opacity=0.8),
                    name=feat,
                    showlegend=False
                ), row=row, col=col)
            else:
                fig.add_trace(go.Scatter3d(
                    x=[0], y=[0], z=[0], text=["No Data"], mode="text",
                    showlegend=False
                ), row=row, col=col)

    fig.update_layout(
        title=dict(
            text="<span style='font-size:32px; font-weight:bold'>Elektrostatik ve Atomik Özellikler</span>",
            x=0.5, xanchor="center", yanchor="top"
        ),
        margin=dict(l=20, r=20, t=90, b=20),
        width=500 * ncols,
        height=500 * nrows,
        paper_bgcolor="white",
        showlegend=False
    )
    fig.update_layout(
        grid=dict(rows=nrows, columns=ncols, pattern="independent"),
        autosize=False,
    )
    for i in range(1, total + 1):
        scene_name = f"scene{i}" if i > 1 else "scene"
        fig.update_layout({scene_name: dict(
            aspectmode="cube",
            xaxis=dict(range=[0, 160], showgrid=True, zeroline=False, showticklabels=False),
            yaxis=dict(range=[0, 160], showgrid=True, zeroline=False, showticklabels=False),
            zaxis=dict(range=[0, 160], showgrid=True, zeroline=False, showticklabels=False),
        )})
    fig.update_traces(showlegend=False)
    fig.show()

def analyze_protein_features(prefix, protein_id, mask_type="surface"):
    h5_path = os.path.join(prefix, f"{protein_id}.h5")
    print(f"\n--- {h5_path} ---")
    with h5py.File(h5_path, "r") as h5f:
        # Summary for each feature and example value
        for feat in h5f["features"]:
            dset = h5f["features"][feat]
            print(f"{feat:20s} shape: {dset.shape} dtype: {dset.dtype}")
            print(f"  min={dset[:].min():.3f} max={dset[:].max():.3f} mean={dset[:].mean():.3f}")
            if dset.ndim == 4:
                print("  Example [0,0,0,0]:", dset[0,0,0,0])
            elif dset.ndim == 3:
                print("  Example [0,0,0]:", dset[0,0,0])
        label = h5f["label/binding_site_calculated"][:]
        label_provided = h5f["label/binding_site_in_dataset"][:]
        print(f"binding_site label shape: {label.shape}, unique values: {np.unique(label)}")
        voxel = (32, 32, 32)
        for feat in h5f["features"]:
            try:
                print(f"{feat:20s}: {h5f['features'][feat][voxel]}")
            except Exception:
                pass

        # Automatic surface & binding/non-binding comparison
        shape_mask = h5f["features/shape"][:]
        dist2surface = h5f["features/dist_to_surface"][:]
        apbs = h5f["features/electrostatic_grid"][:]
        hydrophobicity = h5f["features/hydrophobicity"][:]
        binding = label
        surface_mask = (dist2surface < 2.0).astype(np.uint8)
        shape_mask_bin = (shape_mask > 0).astype(np.uint8)
        intersection = (surface_mask & shape_mask_bin).sum()
        union = (surface_mask | shape_mask_bin).sum()
        iou = intersection / union if union > 0 else np.nan
        dice = (2 * intersection) / (surface_mask.sum() + shape_mask_bin.sum()) if (surface_mask.sum() + shape_mask_bin.sum()) > 0 else np.nan
        print(f"Protein surface: shape vs dist2surface mask IoU={iou:.3f}, Dice={dice:.3f}")
        print(f"  surface voxel count: {surface_mask.sum()} (dist2surface<2A), shape mask voxel count: {shape_mask_bin.sum()}")

        mean_dist2surf_binding = dist2surface[binding==1].mean()
        mean_dist2surf_nonbinding = dist2surface[binding==0].mean()
        print(f"Mean dist2surface (binding): {mean_dist2surf_binding:.2f}")
        print(f"Mean dist2surface (non-binding): {mean_dist2surf_nonbinding:.2f}")
        mean_apbs_binding = apbs[binding==1].mean()
        mean_apbs_nonbinding = apbs[binding==0].mean()
        print(f"Mean APBS (binding): {mean_apbs_binding:.2f}")
        print(f"Mean APBS (non-binding): {mean_apbs_nonbinding:.2f}")
        mean_hydro_binding = hydrophobicity[binding==1].mean()
        mean_hydro_nonbinding = hydrophobicity[binding==0].mean()
        print(f"Mean hydrophobicity (binding): {mean_hydro_binding:.2f}")
        print(f"Mean hydrophobicity (non-binding): {mean_hydro_nonbinding:.2f}")

    print("\n--- Histogram: All features (with mask) ---")
    plot_all_feature_histograms(h5_path, mask_type=mask_type)
    print("\n--- Plotly: Main 3D isosurfaces ---")
    plot_all_features_in_grid(h5_path, atomic_marker_size=3,ncols=6)
    
# USAGE:
#PREFIX = "/Users/tevfik/Sandbox/github/PHD/data/scPDB__cache/box72"
#PREFIX = "/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set-only-used-in-codon-tests/box161"
PREFIX = "/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set_minimal_cache_only_fits_60_resolution_fix/box72"
protein_id = "1a9m"
#protein_id = "1w5y"

analyze_protein_features(PREFIX, protein_id, mask_type="surface")