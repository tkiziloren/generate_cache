import h5py
import numpy as np
import matplotlib.pyplot as plt
import os

def inspect_h5_features(prefix, protein_id):
    path = os.path.join(prefix, f"{protein_id}.h5")
    print(f"--- {path} ---")
    with h5py.File(path, "r") as h5f:
        # Feature özetleri
        for feat in h5f["features"]:
            dset = h5f["features"][feat]
            print(f"{feat:20s} shape: {dset.shape} dtype: {dset.dtype}")
            print(f"  min={dset[:].min():.3f} max={dset[:].max():.3f} mean={dset[:].mean():.3f}")
            if dset.ndim == 4:
                print("  Example [0,0,0,0]:", dset[0,0,0,0])
            elif dset.ndim == 3:
                print("  Example [0,0,0]:", dset[0,0,0])
        # Label
        label = h5f["label/binding_site"][:]
        print(f"binding_site label shape: {label.shape}, unique values: {np.unique(label)}")
        voxel = (32, 32, 32)
        print("APBS:", h5f["features/apbs"][voxel])
        print("Atomic_N:", h5f["features/atomic_N"][voxel])
        print("Dist2Ligand:", h5f["features/dist2ligand"][voxel])
        print("Dist2Surface:", h5f["features/dist2surface"][voxel])
        print("Hydrophobicity:", h5f["features/hydrophobicity"][voxel])

        # ------------------------- MASK ANALİZ -------------------------
        print("\n# --- Yüzey maskesi ve grid alignment otomatik ölçüm ---")
        shape_mask = h5f["features/shape"][:]
        dist2surface = h5f["features/dist2surface"][:]
        apbs = h5f["features/apbs"][:]
        hydrophobicity = h5f["features/hydrophobicity"][:]
        binding = label

        # IoU / Dice Score hesapla (protein yüzeyi için)
        surface_mask = (dist2surface < 2.0).astype(np.uint8)   # 2A yakın voxeller yüzey say
        shape_mask_bin = (shape_mask > 0).astype(np.uint8)
        intersection = (surface_mask & shape_mask_bin).sum()
        union = (surface_mask | shape_mask_bin).sum()
        iou = intersection / union if union > 0 else np.nan
        dice = (2 * intersection) / (surface_mask.sum() + shape_mask_bin.sum()) if (surface_mask.sum() + shape_mask_bin.sum()) > 0 else np.nan
        print(f"Protein yüzeyi: shape vs dist2surface mask IoU={iou:.3f}, Dice={dice:.3f}")
        print(f"  surface voxel sayısı: {surface_mask.sum()} (dist2surface<2A), shape mask voxel sayısı: {shape_mask_bin.sum()}")

        # Binding site vs non-binding site: dist2surface, apbs, hydrophobicity farkları
        print("\n# --- Binding site vs non-binding site karşılaştırma ---")
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

        # ------------------- Histogram: Surface noktaları feature dağılımı -------------------
        try:
            fig, axs = plt.subplots(1, 3, figsize=(12, 3))
            axs[0].hist(apbs[surface_mask==1].flatten(), bins=40, alpha=0.7)
            axs[0].set_title('APBS @ surface')
            axs[1].hist(hydrophobicity[surface_mask==1].flatten(), bins=40, alpha=0.7)
            axs[1].set_title('Hydrophobicity @ surface')
            axs[2].hist(dist2surface.flatten(), bins=40, alpha=0.7)
            axs[2].set_title('Dist2Surface genel')
            plt.tight_layout()
            plt.show()
        except Exception as e:
            print("Histogram çiziminde hata:", e)

    print()
    
import numpy as np
import plotly.graph_objects as go

def plot_3d_isosurface(h5_path, feature, isovalue=None, surface_count=2):
    with h5py.File(h5_path, 'r') as h5f:
        data = h5f['features'][feature][:]
        if isovalue is None:
            isovalue = np.percentile(data, 99)  # otomatik olarak üst %1 eşiği seç
        fig = go.Figure(data=go.Isosurface(
            x=np.arange(data.shape[0]).repeat(data.shape[1]*data.shape[2]),
            y=np.tile(np.arange(data.shape[1]).repeat(data.shape[2]), data.shape[0]),
            z=np.tile(np.arange(data.shape[2]), data.shape[0]*data.shape[1]),
            value=data.flatten(),
            isomin=isovalue, isomax=data.max(),
            surface_count=surface_count,
            colorscale='Viridis',
            opacity=0.6
        ))
        fig.update_layout(title=f"3D Isosurface: {feature}", width=600, height=600)
        fig.show()

# Kullanım örneği:
h5_path = "/Users/tevfik/Sandbox/github/PHD/data/scPDB_minimal_cache/1gzu_1.h5"
plot_3d_isosurface(h5_path, 'apbs')
plot_3d_isosurface(h5_path, 'hydrophobicity')
plot_3d_isosurface(h5_path, 'dist2surface')

# PREFIX = "/Users/tevfik/Sandbox/github/PHD/data/scPDB_minimal_cache"
# for protein_id in ["1gzu_1", "1h0v_1"]:
#     inspect_h5_features(PREFIX, protein_id)
