import argparse
import gzip
import os
import shutil
import subprocess
import tarfile
import urllib.request
import zipfile
from pathlib import Path


DEFAULT_DATA_ROOT = Path("/nfs/production/arl/chembl/tevfik/DEEP_APBS_DATASETS")
DEFAULT_SCPDB_URL = "https://drugdesign.unistra.fr/scPDB/ressources/2016/scPDB.tar.gz"
DEFAULT_P2RANK_REPO = "https://github.com/rdk/p2rank-datasets.git"
DEFAULT_PURESNET_REPO = "https://github.com/jivankandel/PUResNet.git"
DEFAULT_PDBBIND2020_REFINED_URL = ""
PDBBIND_MIN_ARCHIVE_BYTES = 10 * 1024 * 1024


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Stage raw benchmark datasets on Codon/NFS. The script is intentionally "
            "dependency-light so it can run on the datamover partition."
        )
    )
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--force", action="store_true", help="Replace existing prepared directories.")
    parser.add_argument("--skip-scpdb", action="store_true")
    parser.add_argument("--skip-pdbbind", action="store_true")
    parser.add_argument("--skip-puresnet", action="store_true")
    parser.add_argument("--skip-p2rank", action="store_true")
    return parser.parse_args()


def env_url(name, default=""):
    value = os.environ.get(name, "").strip()
    return value or default


def log(message):
    print(message, flush=True)


def ensure_dir(path):
    path.mkdir(parents=True, exist_ok=True)


def local_path_exists(path):
    return path.exists() or path.is_symlink()


def remove_if_force(path, force):
    if local_path_exists(path) and force:
        if path.is_symlink() or path.is_file():
            path.unlink()
        else:
            shutil.rmtree(path)


def download(url, destination):
    if destination.exists() and destination.stat().st_size > 0:
        log(f"[SKIP] archive already exists: {destination}")
        return destination
    if not url:
        raise RuntimeError(f"No URL provided for archive: {destination.name}")
    ensure_dir(destination.parent)
    tmp = destination.with_suffix(destination.suffix + ".part")
    log(f"[GET] {url}")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(destination)
    return destination


def find_existing_archive(archive_dir, names):
    for name in names:
        candidate = archive_dir / name
        if candidate.exists() and candidate.stat().st_size > 0:
            return candidate
    return None


def validate_tar_archive(archive, label, min_size_bytes=0):
    size = archive.stat().st_size
    if min_size_bytes and size < min_size_bytes:
        raise RuntimeError(
            f"{label} archive is too small ({size} bytes): {archive}. "
            "This usually means the download returned an HTML/login/error page, "
            "not the real dataset archive."
        )
    if not tarfile.is_tarfile(archive):
        raise RuntimeError(
            f"{label} archive is not a valid tar file: {archive}. "
            "Place the official tar.gz archive under the archives directory and rerun."
        )


def safe_extract_tar(tf, destination):
    destination_resolved = destination.resolve()
    for member in tf.getmembers():
        member_path = (destination / member.name).resolve()
        try:
            member_path.relative_to(destination_resolved)
        except ValueError as exc:
            raise RuntimeError(f"Unsafe path in archive: {member.name}")
    tf.extractall(destination)


def safe_extract_zip(zf, destination):
    destination_resolved = destination.resolve()
    for member in zf.infolist():
        member_path = (destination / member.filename).resolve()
        try:
            member_path.relative_to(destination_resolved)
        except ValueError as exc:
            raise RuntimeError(f"Unsafe path in archive: {member.filename}")
    zf.extractall(destination)


def extract_archive(archive, destination, force=False):
    remove_if_force(destination, force)
    if destination.exists() and any(destination.iterdir()):
        log(f"[SKIP] already extracted: {destination}")
        return destination
    ensure_dir(destination)
    log(f"[EXTRACT] {archive} -> {destination}")
    name = archive.name.lower()
    if name.endswith(".zip"):
        with zipfile.ZipFile(archive) as zf:
            safe_extract_zip(zf, destination)
    elif name.endswith((".tar.gz", ".tgz", ".tar")):
        try:
            with tarfile.open(archive) as tf:
                safe_extract_tar(tf, destination)
        except tarfile.TarError as exc:
            raise RuntimeError(
                f"Could not extract {archive}. If this is a PDBBind archive, "
                "the download may require licensed/manual access. Place the "
                "official tar.gz under the archives directory and rerun."
            ) from exc
    elif name.endswith(".gz"):
        out = destination / archive.with_suffix("").name
        with gzip.open(archive, "rb") as src, out.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        raise RuntimeError(f"Unsupported archive type: {archive}")
    return destination


def move_contents(src, dst, force=False):
    remove_if_force(dst, force)
    ensure_dir(dst)
    for child in src.iterdir():
        target = dst / child.name
        if target.exists():
            continue
        shutil.move(str(child), str(target))


def find_scpdb_payload(extracted_root):
    candidates = [extracted_root]
    candidates.extend(path for path in extracted_root.rglob("*") if path.is_dir())
    for candidate in candidates:
        case_dirs = list(candidate.glob("*"))
        hits = [
            path
            for path in case_dirs[:500]
            if path.is_dir() and (path / "protein.mol2").exists() and (path / "ligand.mol2").exists()
        ]
        if hits:
            return candidate
    raise RuntimeError(f"Could not find scPDB case directories under {extracted_root}")


def find_pdbbind_refined_payload(extracted_root):
    direct = extracted_root / "refined-set"
    if direct.exists():
        return direct
    candidates = [path for path in extracted_root.rglob("*") if path.is_dir()]
    for candidate in candidates:
        case_dirs = [path for path in candidate.iterdir() if path.is_dir()]
        checked = case_dirs[:500]
        hits = []
        for case_dir in checked:
            case = case_dir.name
            if (case_dir / f"{case}_protein.pdb").exists() and (
                (case_dir / f"{case}_ligand.mol2").exists()
                or (case_dir / f"{case}_ligand.sdf").exists()
            ):
                hits.append(case_dir)
        if len(hits) >= 5:
            return candidate
    raise RuntimeError(f"Could not find PDBBind refined-set case directories under {extracted_root}")


def clone_or_pull(repo_url, destination):
    if destination.exists() and (destination / ".git").exists():
        log(f"[PULL] {destination}")
        subprocess.run(["git", "-C", str(destination), "pull", "--ff-only"], check=True)
        return
    if destination.exists():
        raise RuntimeError(f"Destination exists but is not a git repository: {destination}")
    ensure_dir(destination.parent)
    log(f"[CLONE] {repo_url} -> {destination}")
    subprocess.run(["git", "clone", repo_url, str(destination)], check=True)


def prepare_scpdb(data_root, archives, staging, force):
    url = env_url("SCPDB_URL", DEFAULT_SCPDB_URL)
    archive = find_existing_archive(archives, ["scPDB.tar.gz", "scPDB.tgz", "scPDB.tar"])
    if archive is None:
        archive = download(url, archives / "scPDB.tar.gz")
    extracted = extract_archive(archive, staging / "scpdb_extract", force=force)
    payload = find_scpdb_payload(extracted)
    final = data_root / "datasets" / "scPDB"
    move_contents(payload, final, force=force)
    log(f"[OK] scPDB -> {final}")


def prepare_pdbbind(data_root, archives, staging, force):
    refined_url = env_url("PDBBIND_REFINED_URL", "")
    refined_2020_url = env_url("PDBBIND2020_REFINED_URL", DEFAULT_PDBBIND2020_REFINED_URL)

    archive = find_existing_archive(
        archives,
        [
            "pdbbind_refined.tar.gz",
            "pdbbind_refined.tgz",
            "pdbbind_v2020_refined.tar.gz",
            "PDBbind_v2020_refined.tar.gz",
        ],
    )
    if archive is None:
        if not refined_2020_url:
            raise RuntimeError(
                "PDBBind refined-set archive was not found under the archives directory. "
                "PDBBind usually requires licensed/manual access, so this script does not "
                "guess a public URL. Place pdbbind_v2020_refined.tar.gz under archives/ "
                "or provide PDBBIND2020_REFINED_URL explicitly."
            )
        archive = download(refined_2020_url, archives / "pdbbind_v2020_refined.tar.gz")

    validate_tar_archive(archive, "PDBBind refined-set", PDBBIND_MIN_ARCHIVE_BYTES)
    extracted = extract_archive(archive, staging / "pdbbind2020_refined_extract", force=force)
    payload = find_pdbbind_refined_payload(extracted)
    final_2020 = data_root / "datasets" / "pdbbind2020" / "refined-set"
    move_contents(payload, final_2020, force=force)
    log(f"[OK] PDBBind 2020 refined-set -> {final_2020}")

    legacy_final = data_root / "datasets" / "pdbbind" / "refined-set"
    if refined_url:
        legacy_archive = download(refined_url, archives / "pdbbind_refined.tar.gz")
        legacy_extracted = extract_archive(legacy_archive, staging / "pdbbind_refined_extract", force=force)
        legacy_payload = find_pdbbind_refined_payload(legacy_extracted)
        move_contents(legacy_payload, legacy_final, force=force)
    elif not local_path_exists(legacy_final):
        ensure_dir(legacy_final.parent)
        rel_target = os.path.relpath(final_2020, legacy_final.parent)
        legacy_final.symlink_to(rel_target)
        log(f"[LINK] {legacy_final} -> {rel_target}")
    log(f"[OK] PDBBind legacy alias -> {legacy_final}")


def prepare_puresnet(data_root, force):
    repo = data_root / "sources" / "PUResNet"
    clone_or_pull(env_url("PURESNET_REPO", DEFAULT_PURESNET_REPO), repo)
    target = data_root / "datasets" / "external_benchmarks" / "puresnet"
    ensure_dir(target)
    for zip_name, subdir in (
        ("coach.zip", "coach"),
        ("BU48.zip", "BU48"),
        ("scpdb_subset.zip", "scpdb_subset"),
    ):
        archive = repo / zip_name
        if archive.exists():
            extract_archive(archive, target / subdir, force=force)
        else:
            log(f"[WARN] Missing PUResNet archive in repository: {archive}")
    log(f"[OK] PUResNet benchmark files -> {target}")


def prepare_p2rank(data_root):
    destination = data_root / "datasets" / "external_benchmarks" / "p2rank-datasets"
    clone_or_pull(env_url("P2RANK_REPO", DEFAULT_P2RANK_REPO), destination)
    log(f"[OK] P2Rank datasets -> {destination}")


def write_layout(data_root):
    layout = data_root / "DATASET_LAYOUT.md"
    text = f"""# DEEP_APBS_DATASETS

Root:

```text
{data_root}
```

Raw archives:

```text
{data_root / "archives"}
```

Prepared dataset roots:

```text
{data_root / "datasets" / "scPDB"}
{data_root / "datasets" / "pdbbind" / "refined-set"}
{data_root / "datasets" / "pdbbind2020" / "refined-set"}
{data_root / "datasets" / "external_benchmarks" / "puresnet"}
{data_root / "datasets" / "external_benchmarks" / "p2rank-datasets"}
```

Work11 cache output root:

```text
{data_root / "cache" / "work11_cache_gridfix_v1"}
```

Notes:

- `datasets/pdbbind/refined-set` is a legacy alias. If no separate legacy PDBBind
  URL is provided, it points to `datasets/pdbbind2020/refined-set`.
- PDBBind may require licensed/manual access. If direct download fails, place the
  archive under `archives/` and rerun the datamover job.
"""
    layout.write_text(text, encoding="utf-8")
    log(f"[OK] wrote {layout}")


def main():
    args = parse_args()
    data_root = args.data_root
    archives = data_root / "archives"
    staging = data_root / "staging"
    ensure_dir(archives)
    ensure_dir(staging)
    ensure_dir(data_root / "datasets")
    ensure_dir(data_root / "cache")

    if not args.skip_scpdb:
        prepare_scpdb(data_root, archives, staging, args.force)
    if not args.skip_pdbbind:
        prepare_pdbbind(data_root, archives, staging, args.force)
    if not args.skip_puresnet:
        prepare_puresnet(data_root, args.force)
    if not args.skip_p2rank:
        prepare_p2rank(data_root)
    write_layout(data_root)
    log("[DONE] dataset staging finished")


if __name__ == "__main__":
    main()
