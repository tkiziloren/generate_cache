import h5py
with h5py.File("/Users/tevfik/Sandbox/github/PHD/data/scPDB_minimal_cache/1a2b_1.h5") as h5f:
    grid = h5f["features/atomic_N"][:]  # (64, 64, 64) veya (65,65,65)
    print(grid[9, 41, 31])