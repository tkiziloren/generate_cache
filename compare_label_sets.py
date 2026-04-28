import h5py
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns  # for nicer heatmap

def compare_label_sets(h5_path):
    """
    İki farklı binding site etiketi için:
    - Histogramları üst üste çiz
    - Confusion matrix (TP, TN, FP, FN), IoU ve Dice skorlarını hesapla
    - Ortadan alınan bir dilim (slice) üzerinde farkları görselleştir
    """
    with h5py.File(h5_path, "r") as f:
        lbl_calc = f["label/binding_site_calculated"][:]
        lbl_data = f["label/binding_site_in_dataset"][:]
    
    # Flatten
    flat_calc = lbl_calc.flatten()
    flat_data = lbl_data.flatten()
    
    # 1) Overlay histogram
    plt.figure(figsize=(6,4))
    plt.hist(flat_calc, bins=[-0.5,0.5,1.5], alpha=0.6, label="Calculated")
    plt.hist(flat_data, bins=[-0.5,0.5,1.5], alpha=0.6, label="Dataset")
    plt.xticks([0,1])
    plt.xlabel("Label Değeri")
    plt.ylabel("Vokseller")
    plt.legend()
    plt.title("Overlay Histogram: Calculated vs Dataset")
    plt.tight_layout()
    plt.show()

    # 2) Confusion matrix
    TP = np.sum((flat_calc==1) & (flat_data==1))
    TN = np.sum((flat_calc==0) & (flat_data==0))
    FP = np.sum((flat_calc==0) & (flat_data==1))
    FN = np.sum((flat_calc==1) & (flat_data==0))
    
    # IoU ve Dice
    iou = TP / (TP + FP + FN) if (TP + FP + FN) > 0 else 0
    dice = 2*TP / (2*TP + FP + FN) if (2*TP + FP + FN) > 0 else 0
    
    print(f"TP: {TP}, TN: {TN}, FP: {FP}, FN: {FN}")
    print(f"IoU: {iou:.4f}, Dice: {dice:.4f}")
    
    # Confusion matrix heatmap
    cm = np.array([[TN, FP],
                   [FN, TP]])
    plt.figure(figsize=(4,4))
    sns.heatmap(cm, annot=True, fmt="d", cbar=False,
                xticklabels=["Calc=0","Calc=1"],
                yticklabels=["Data=0","Data=1"])
    plt.xlabel("Calculated")
    plt.ylabel("Dataset")
    plt.title("Confusion Matrix")
    plt.tight_layout()
    plt.show()
    
    # 3) Voxel-level fark haritası: orta dilim
    shape = lbl_calc.shape
    mid_z = shape[2]//2
    slice_calc = lbl_calc[:,:,mid_z]
    slice_data = lbl_data[:,:,mid_z]
    diff_slice = slice_calc - slice_data  # +1: calc only, -1: data only, 0: same
    
    plt.figure(figsize=(6,6))
    plt.imshow(diff_slice, cmap="bwr", vmin=-1, vmax=1)
    plt.colorbar(ticks=[-1,0,1], label="Fark (Calc - Data)")
    plt.title(f"Difference Map Slice at z={mid_z}")
    plt.tight_layout()
    plt.show()


import h5py
import numpy as np
import matplotlib.pyplot as plt

def analyze_and_plot_labels(h5_path):
    # 1) Veriyi oku
    with h5py.File(h5_path, "r") as f:
        lbl_calc = f["label/binding_site_calculated"][:].flatten()
        lbl_data = f["label/binding_site_in_dataset"][:].flatten()

    # 2) Sayısal metrikler
    TP = np.sum((lbl_calc==1) & (lbl_data==1))
    TN = np.sum((lbl_calc==0) & (lbl_data==0))
    FP = np.sum((lbl_calc==0) & (lbl_data==1))
    FN = np.sum((lbl_calc==1) & (lbl_data==0))
    iou  = TP / (TP + FP + FN) if (TP + FP + FN)>0 else 0
    dice = 2*TP / (2*TP + FP + FN) if (2*TP + FP + FN)>0 else 0
    print(f"TP={TP}, TN={TN}, FP={FP}, FN={FN}")
    print(f"IoU={iou:.4f},  Dice={dice:.4f}\n")

    # 3) Normalize edilmiş ve log­ölçekli overlay histogram
    bins = [-0.5, 0.5, 1.5]
    plt.figure(figsize=(6,4))
    plt.hist(lbl_calc, bins=bins, alpha=0.6, label="Calculated",
             density=True, log=True)
    plt.hist(lbl_data, bins=bins, alpha=0.6, label="Dataset",
             density=True, log=True)
    plt.xticks([0,1])
    plt.xlabel("Label Değeri")
    plt.ylabel("Normalize edilip log ölçek")
    plt.legend()
    plt.title("Overlay (norm. & log): Calc vs Data")
    plt.tight_layout()
    plt.show()

    # 4) Yalnızca pozitif (1) etiketleri karşılaştır
    count_calc_1 = np.sum(lbl_calc==1)
    count_data_1 = np.sum(lbl_data==1)
    plt.figure(figsize=(5,4))
    plt.bar([0,1], [count_calc_1, count_data_1], width=0.4, 
            label=["Calculated","Dataset"], color=["C0","C1"])
    plt.xticks([0,1], ["Calculated","Dataset"])
    plt.ylabel("1 etiketli voksel sayısı")
    plt.title("Pozitif Voksel Sayıları")
    plt.tight_layout()
    plt.show()

    # 5) Orta dilimde fark haritası
    # (1: sadece calc, -1: sadece data, 0: ikisinde de aynı)
    # diyelim 3D shape = (X,Y,Z)
    with h5py.File(h5_path, "r") as f:
        shape3d = f["label/binding_site_calculated"].shape
    z0 = shape3d[2]//2
    slc = (lbl_calc.reshape(shape3d)[:,:,z0] - 
           lbl_data.reshape(shape3d)[:,:,z0])
    plt.figure(figsize=(5,5))
    plt.imshow(slc, cmap="bwr", vmin=-1, vmax=1)
    plt.colorbar(ticks=[-1,0,1], label="Calc–Data farkı")
    plt.title(f"Orta dilim (z={z0}) fark haritası")
    plt.tight_layout()
    plt.show()
# Örnek kullanım:
analyze_and_plot_labels("/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set_minimal_cache_only_fits_60/box72/1a1e.h5")
