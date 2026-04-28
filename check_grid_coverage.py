import h5py
import numpy as np
import prody

# Test için bir h5 dosyası
h5_path = "/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set_minimal/box161/1a1e.h5"
protein_pdb = "/Users/tevfik/Sandbox/github/PHD/data/pdbbind/refined-set_minimal/1a1e/1a1e_protein.pdb"

with h5py.File(h5_path, "r") as h5f:
    # APBS grid
    apbs = h5f['features/electrostatic_grid'][:]
    
    # Atomic C
    atomic_c = h5f['features/atomic_C'][:]
    
    # Metadata
    center = np.array(h5f.attrs['center'])
    box_size = h5f.attrs['box_size']
    resolution = h5f.attrs['resolution']
    
    print(f"Box size: {box_size}")
    print(f"Resolution: {resolution}")
    print(f"Grid center: {center}")
    print(f"\nAPBS shape: {apbs.shape}")
    print(f"Atomic_C shape: {atomic_c.shape}")
    
    # APBS'de sıfırdan farklı değerler
    apbs_nonzero = np.where(apbs != 0)
    if len(apbs_nonzero[0]) > 0:
        apbs_min = apbs_nonzero[0].min(), apbs_nonzero[1].min(), apbs_nonzero[2].min()
        apbs_max = apbs_nonzero[0].max(), apbs_nonzero[1].max(), apbs_nonzero[2].max()
        print(f"\nAPBS non-zero voxel range:")
        print(f"  Min indices: {apbs_min}")
        print(f"  Max indices: {apbs_max}")
        print(f"  Span: {apbs_max[0]-apbs_min[0]+1} x {apbs_max[1]-apbs_min[1]+1} x {apbs_max[2]-apbs_min[2]+1}")
    else:
        print("\nNo non-zero APBS values!")
    
    # Atomic C'de sıfırdan farklı değerler
    atomic_nonzero = np.where(atomic_c > 0)
    if len(atomic_nonzero[0]) > 0:
        atomic_min = atomic_nonzero[0].min(), atomic_nonzero[1].min(), atomic_nonzero[2].min()
        atomic_max = atomic_nonzero[0].max(), atomic_nonzero[1].max(), atomic_nonzero[2].max()
        print(f"\nAtomic_C non-zero voxel range:")
        print(f"  Min indices: {atomic_min}")
        print(f"  Max indices: {atomic_max}")
        print(f"  Span: {atomic_max[0]-atomic_min[0]+1} x {atomic_max[1]-atomic_min[1]+1} x {atomic_max[2]-atomic_min[2]+1}")
    else:
        print("\nNo atomic_C data!")
    
    # Fiziksel boyutlar (Angstrom cinsinden)
    print(f"\n=== FİZİKSEL BOYUTLAR (Angstrom) ===")
    grid_start = center - (box_size // 2) * resolution
    grid_end = center + (box_size // 2) * resolution
    print(f"Grid başlangıç: {grid_start}")
    print(f"Grid bitiş: {grid_end}")
    print(f"Toplam fiziksel alan: {grid_end - grid_start} Å")
    print(f"Beklenen: {(box_size - 1) * resolution} Å (tfbio formülü)")
    
    # APBS değer dağılımı
    print(f"\n=== APBS DEĞER DAĞILIMI ===")
    print(f"APBS min: {apbs.min():.2f}")
    print(f"APBS max: {apbs.max():.2f}")
    print(f"APBS mean: {apbs.mean():.2f}")
    print(f"APBS std: {apbs.std():.2f}")
    
    # -15 ile +15 arasında kaç voxel var?
    in_range = ((apbs >= -15) & (apbs <= 15))
    print(f"\nAPBS değerleri -15 ile +15 arasında: {in_range.sum()} / {apbs.size} voxel")
    
    # Farklı percentile'lar
    percentiles = [1, 5, 10, 25, 50, 75, 90, 95, 99]
    print(f"\nAPBS Percentiles:")
    for p in percentiles:
        val = np.percentile(apbs, p)
        print(f"  {p:2d}%: {val:8.2f}")
    
    # Shape mask ile karşılaştır
    shape = h5f['features/shape'][:]
    protein_region = (shape > 0.5)
    apbs_in_protein = apbs[protein_region]
    print(f"\n=== APBS (Sadece Protein Bölgesinde) ===")
    print(f"Protein voxel sayısı: {protein_region.sum()}")
    if apbs_in_protein.size > 0:
        print(f"APBS (protein) min: {apbs_in_protein.min():.2f}")
        print(f"APBS (protein) max: {apbs_in_protein.max():.2f}")
        print(f"APBS (protein) mean: {apbs_in_protein.mean():.2f}")
    
    # Shape coverage
    shape_nonzero = np.where(shape > 0.5)
    if len(shape_nonzero[0]) > 0:
        shape_min = shape_nonzero[0].min(), shape_nonzero[1].min(), shape_nonzero[2].min()
        shape_max = shape_nonzero[0].max(), shape_nonzero[1].max(), shape_nonzero[2].max()
        print(f"\n=== SHAPE (protein surface) ===")
        print(f"  Min indices: {shape_min}")
        print(f"  Max indices: {shape_max}")
        print(f"  Span: {shape_max[0]-shape_min[0]+1} x {shape_max[1]-shape_min[1]+1} x {shape_max[2]-shape_min[2]+1}")
        print(f"  Total voxels: {protein_region.sum()}")
        
    print(f"\n=== KARŞILAŞTIRMA ===")
    print(f"Shape span:    {shape_max[0]-shape_min[0]+1} x {shape_max[1]-shape_min[1]+1} x {shape_max[2]-shape_min[2]+1}")
    print(f"Atomic_C span: {atomic_max[0]-atomic_min[0]+1} x {atomic_max[1]-atomic_min[1]+1} x {atomic_max[2]-atomic_min[2]+1}")
    print(f"\nShape ve Atomic_C AYNI OLMALI ama farklılar!")
    print(f"Shape büyüklüğü: ~{np.mean([shape_max[i]-shape_min[i]+1 for i in range(3)]):.1f} voxel")
    print(f"Atomic_C büyüklüğü: ~{np.mean([atomic_max[i]-atomic_min[i]+1 for i in range(3)]):.1f} voxel")
    
    # Gerçek protein koordinatlarını kontrol et
    print(f"\n=== GERÇEK PROTEİN KOORDİNATLARI ===")
    protein = prody.parsePDB(protein_pdb)
    protein_coords = protein.getCoords()
    print(f"Protein atom sayısı: {len(protein_coords)}")
    print(f"Protein coord min: {protein_coords.min(axis=0)}")
    print(f"Protein coord max: {protein_coords.max(axis=0)}")
    print(f"Protein coord center: {protein_coords.mean(axis=0)}")
    print(f"Grid center (h5 metadata): {center}")
    print(f"\nFark: {protein_coords.mean(axis=0) - center}")
    
    # Shape mask nasıl hesaplanmış?
    print(f"\n=== SHAPE MASK HESAPLAMA TESTİ ===")
    grid_start = center - (box_size // 2) * resolution
    print(f"Grid start: {grid_start}")
    
    # İlk carbon atomunun indexini hesapla (manual)
    carbon_atoms = protein.select('element C')
    if carbon_atoms:
        first_c = carbon_atoms.getCoords()[0]
        expected_idx = np.round((first_c - grid_start) / resolution).astype(int)
        print(f"İlk C atomu gerçek coord: {first_c}")
        print(f"Beklenen grid index: {expected_idx}")
        print(f"Atomic_C grid'de bu noktada değer var mı? {atomic_c[tuple(expected_idx)] if np.all((expected_idx >= 0) & (expected_idx < box_size)) else 'OUT OF BOUNDS!'}")
