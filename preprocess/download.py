"""Download helpers for public data sources used by PEREGRINE."""

import tarfile
import gzip
import shutil
from pathlib import Path
from urllib.request import urlretrieve, urlopen
from urllib.error import URLError, HTTPError


GTEX_V7_EQTL_URL = (
    "https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/"
    "GTEx_Analysis_v7_eQTL.tar.gz"
)

PSYCHIC_BASE = "https://www.cs.huji.ac.il/~tommy/PSYCHIC/"

# Human hierarchical TAD enhancer files at p < 1e-4 used by original workflow
PSYCHIC_FILES = [
    "hES.enh_1e-4.bed",
    "hIMR90.enh_1e-4.bed",
    "lib25_GM12878.enh_1e-4.bed",
    "lib25_HMEC.enh_1e-4.bed",
    "lib25_HUVEC.enh_1e-4.bed",
    "lib25_IMR90.enh_1e-4.bed",
    "lib25_K562.enh_1e-4.bed",
    "lib25_KBM7.enh_1e-4.bed",
    "lib25_NHEK.enh_1e-4.bed",
    "mCO.enh_1e-4.bed",
    "mES.enh_1e-4.bed",
    "pombo_ESC_HindIII.enh_1e-4.bed",
    "pombo_ESC_NcoI.enh_1e-4.bed",
    "pombo_Neuron_HindIII.enh_1e-4.bed",
    "pombo_Neuron_NcoI.enh_1e-4.bed",
    "pombo_NPC_HindIII.enh_1e-4.bed",
    "pombo_NPC_NcoI.enh_1e-4.bed",
]

ENCODE_DOMAIN_FILES = {
    "ENCFF274VJU.bed": "https://www.encodeproject.org/files/ENCFF274VJU/@@download/ENCFF274VJU.bed.gz",
    "ENCFF588KUZ.bed": "https://www.encodeproject.org/files/ENCFF588KUZ/@@download/ENCFF588KUZ.bed.gz",
}


def _download(url: str, dest: Path, force: bool = False):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not force:
        print(f"Already present: {dest}")
        return dest
    print(f"Downloading {url} -> {dest}")
    try:
        urlretrieve(url, dest)
    except (URLError, HTTPError) as exc:
        raise RuntimeError(f"Failed to download {url}: {exc}") from exc
    return dest


def download_gtex_v7(target_dir: Path, force: bool = False):
    """Download and extract GTEx v7 significant variant-gene pairs."""
    target_dir = Path(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)
    archive = target_dir / "GTEx_Analysis_v7_eQTL.tar.gz"
    _download(GTEX_V7_EQTL_URL, archive, force=force)

    extract_dir = target_dir / "GTEx_Analysis_v7_eQTL"
    extract_dir.mkdir(exist_ok=True)

    print("Extracting significant variant-gene pair files...")
    with tarfile.open(archive, "r:gz") as tar:
        members = [
            m for m in tar.getmembers()
            if m.name.endswith("signif_variant_gene_pairs.txt.gz")
        ]
        if not members:
            raise RuntimeError("No *signif_variant_gene_pairs.txt.gz files found in GTEx archive.")
        for member in members:
            # Flatten into extract_dir
            member.name = Path(member.name).name
            tar.extract(member, path=extract_dir)
    print(f"GTEx files ready in {extract_dir}")
    return extract_dir


def download_psychic(target_dir: Path, force: bool = False):
    """Download PSYCHIC hierarchical TAD enhancer BED files (p < 1e-4)."""
    target_dir = Path(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)
    downloaded = []
    for fname in PSYCHIC_FILES:
        url = PSYCHIC_BASE + fname
        dest = target_dir / fname
        try:
            _download(url, dest, force=force)
            downloaded.append(dest)
        except RuntimeError as exc:
            print(f"Warning: {exc}")
    if not downloaded:
        raise RuntimeError("No PSYCHIC files could be downloaded.")
    print(f"Downloaded {len(downloaded)} PSYCHIC BED files to {target_dir}")
    return downloaded


def download_encode_domains(target_dir: Path, force: bool = False):
    """Download ENCODE nested contact domain / TAD BED files used by the original pipeline."""
    target_dir = Path(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for out_name, url in ENCODE_DOMAIN_FILES.items():
        gz_dest = target_dir / (out_name + ".gz")
        final_dest = target_dir / out_name
        if final_dest.exists() and not force:
            print(f"Already present: {final_dest}")
            results.append(final_dest)
            continue
        _download(url, gz_dest, force=force)
        with gzip.open(gz_dest, "rb") as src, open(final_dest, "wb") as dst:
            shutil.copyfileobj(src, dst)
        gz_dest.unlink(missing_ok=True)
        results.append(final_dest)
    print(f"ENCODE domain files ready in {target_dir}")
    return results
