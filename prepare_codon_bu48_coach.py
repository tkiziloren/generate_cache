import argparse
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path


DEFAULT_DATA_ROOT = Path("/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS")
DEFAULT_P2RANK_REPO = "https://github.com/rdk/p2rank-datasets.git"
DEFAULT_PURESNET_REPO = "https://github.com/jivankandel/PUResNet.git"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Stage only BU48 and COACH420 external benchmarks on Codon/NFS."
    )
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--repo-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--skip-prepared", action="store_true")
    parser.add_argument("--p2rank-repo", default=DEFAULT_P2RANK_REPO)
    parser.add_argument("--puresnet-repo", default=DEFAULT_PURESNET_REPO)
    return parser.parse_args()


def log(message):
    print(message, flush=True)


def ensure_dir(path):
    path.mkdir(parents=True, exist_ok=True)


def remove_if_force(path, force):
    if path.exists() and force:
        if path.is_symlink() or path.is_file():
            path.unlink()
        else:
            shutil.rmtree(path)


def run(cmd):
    log("[RUN] " + " ".join(str(part) for part in cmd))
    subprocess.run([str(part) for part in cmd], check=True)


def clone_or_pull(repo_url, destination):
    if destination.exists() and (destination / ".git").exists():
        log(f"[PULL] {destination}")
        run(["git", "-C", destination, "pull", "--ff-only"])
        return
    if destination.exists():
        raise RuntimeError(f"Destination exists but is not a git repository: {destination}")
    ensure_dir(destination.parent)
    log(f"[CLONE] {repo_url} -> {destination}")
    run(["git", "clone", repo_url, destination])


def sparse_clone_or_pull(repo_url, destination, sparse_paths, force=False):
    remove_if_force(destination, force)
    if destination.exists() and (destination / ".git").exists():
        log(f"[SPARSE PULL] {destination}")
        run(["git", "-C", destination, "sparse-checkout", "init", "--cone"])
        run(["git", "-C", destination, "sparse-checkout", "set", *sparse_paths])
        run(["git", "-C", destination, "pull", "--ff-only"])
        return
    if destination.exists():
        raise RuntimeError(f"Destination exists but is not a git repository: {destination}")
    ensure_dir(destination.parent)
    log(f"[SPARSE CLONE] {repo_url} -> {destination}")
    run(["git", "clone", "--depth", "1", "--filter=blob:none", "--sparse", repo_url, destination])
    run(["git", "-C", destination, "sparse-checkout", "set", *sparse_paths])


def safe_extract_zip(archive, destination, force=False):
    remove_if_force(destination, force)
    if destination.exists() and any(destination.iterdir()):
        log(f"[SKIP] already extracted: {destination}")
        return
    ensure_dir(destination)
    destination_resolved = destination.resolve()
    log(f"[EXTRACT] {archive} -> {destination}")
    with zipfile.ZipFile(archive) as zf:
        for member in zf.infolist():
            member_path = (destination / member.filename).resolve()
            try:
                member_path.relative_to(destination_resolved)
            except ValueError as exc:
                raise RuntimeError(f"Unsafe path in archive: {member.filename}") from exc
        zf.extractall(destination)


def require_file(path, label):
    if not path.exists() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty {label}: {path}")


def prepare_puresnet_sources(data_root, repo_url, force):
    sources_root = data_root / "sources"
    repo = sources_root / "PUResNet"
    clone_or_pull(repo_url, repo)

    puresnet_root = data_root / "datasets" / "external_benchmarks" / "puresnet"
    ensure_dir(puresnet_root)
    coach_zip = repo / "coach.zip"
    bu48_zip = repo / "BU48.zip"
    require_file(coach_zip, "PUResNet coach.zip")
    require_file(bu48_zip, "PUResNet BU48.zip")
    safe_extract_zip(coach_zip, puresnet_root / "coach", force=force)
    safe_extract_zip(bu48_zip, puresnet_root / "BU48", force=force)
    log(f"[OK] PUResNet BU48/COACH -> {puresnet_root}")
    return puresnet_root


def prepare_p2rank_subset(data_root, repo_url, force):
    p2rank_root = data_root / "datasets" / "external_benchmarks" / "p2rank-datasets"
    sparse_clone_or_pull(
        repo_url,
        p2rank_root,
        sparse_paths=["coach420", "joined/bu48", "_lists"],
        force=force,
    )
    log(f"[OK] P2Rank BU48/COACH subset -> {p2rank_root}")
    return p2rank_root


def prepare_matched_benchmark(repo_dir, data_root, puresnet_root, p2rank_root, force):
    output_root = data_root / "datasets" / "external_benchmarks" / "puresnet_prepared"
    remove_if_force(output_root, force)
    script = repo_dir / "prepare_puresnet_benchmarks.py"
    obabel = repo_dir / ".conda" / "bin" / "obabel"
    require_file(script, "prepare_puresnet_benchmarks.py")
    require_file(obabel, "Open Babel executable")
    cmd = [
        sys.executable,
        script,
        "--puresnet-root",
        puresnet_root,
        "--p2rank-root",
        p2rank_root,
        "--output-root",
        output_root,
        "--obabel-bin",
        obabel,
    ]
    if force:
        cmd.append("--overwrite")
    run(cmd)
    log(f"[OK] prepared benchmark subset -> {output_root}")


def write_layout(data_root):
    layout = data_root / "BU48_COACH_LAYOUT.md"
    text = f"""# BU48 and COACH420 dataset layout

Prepared under:

```text
{data_root / "datasets" / "external_benchmarks"}
```

Expected source directories:

```text
{data_root / "datasets" / "external_benchmarks" / "puresnet" / "coach"}
{data_root / "datasets" / "external_benchmarks" / "puresnet" / "BU48"}
{data_root / "datasets" / "external_benchmarks" / "p2rank-datasets" / "coach420"}
{data_root / "datasets" / "external_benchmarks" / "p2rank-datasets" / "joined" / "bu48"}
```

Prepared matched benchmark directory:

```text
{data_root / "datasets" / "external_benchmarks" / "puresnet_prepared"}
```
"""
    layout.write_text(text, encoding="utf-8")
    log(f"[OK] wrote {layout}")


def main():
    args = parse_args()
    ensure_dir(args.data_root / "datasets" / "external_benchmarks")
    puresnet_root = prepare_puresnet_sources(args.data_root, args.puresnet_repo, args.force)
    p2rank_root = prepare_p2rank_subset(args.data_root, args.p2rank_repo, args.force)
    if not args.skip_prepared:
        prepare_matched_benchmark(args.repo_dir, args.data_root, puresnet_root, p2rank_root, args.force)
    write_layout(args.data_root)
    log("[DONE] BU48/COACH staging finished")


if __name__ == "__main__":
    main()
