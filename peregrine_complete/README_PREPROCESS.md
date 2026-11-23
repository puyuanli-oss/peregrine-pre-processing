# PEREGRINE preprocessing (extended)

Run `python -m unittest -v tests_preprocessing` for local verification.

**Complete preprocessing:** `./run_preprocess.sh all`

The `all` command downloads GTEx v7, PSYCHIC, ENCODE, and a pinned Ensembl GRCh37 release-87 GTF, then generates:
- `preprocessed/promoters_600bp.bed` (strand-aware 600-bp promoters)
- `eQTL/exons_genes.txt` (exon annotation)
- `tad/TSSgenesbed` (TSS annotation)
- `preprocessed/eqtl/*.bed` (one eQTL intermediate per GTEx tissue)
- `preprocessed/tad/*.bed` (cleaned ENCODE domains)
- cleaned PSYCHIC BEDs via the existing cleanup helper

**Offline / cached mode:** `./run_preprocess.sh all --offline --gtf /path/to/GRCh37.gtf.gz`
Place GTEx files in `eQTL/GTEx_Analysis_v7_eQTL/`, PSYCHIC files in `heirarchicalTAD/data/`, and ENCODE files in `tad/` beforehand. `--steps gtex,encode` selects only those downloaded-data transformations; annotation processing always runs. `--use-panther` restricts promoters to PANTHER-listed genes.

**Important scientific/operational caveats:** GRCh37 is required to match the original hg19-based data; other GTF builds are not automatically detected. GTF gene features supply a single gene-level TSS, which is **not identical** to a full transcript-level BioMart TSS export. The original `tad/TSSgenesbed` is replaced when `all` runs; back it up if you need the historical version. The original ChIA-PET and other linking scripts are retained but are **not** run by `all`. The biocurated `CREbedDBenhancers_10092018` is absent from the supplied archive and must be supplied separately for downstream linking and the original `check` command. Large external downloads and full downstream scientific output have not been tested in this environment. The existing PSYCHIC downloader tolerates partial availability; inspect tissue coverage before interpreting results.

The original README and the rest of the original scripts are preserved. The new `preprocess/auto.py` and `tests_preprocessing.py` provide the extensions.

---

# PEREGRINE Pre-processing Pipeline

This directory adds an automated pre-processing layer that produces the intermediate files expected by the existing PEREGRINE linking scripts (eQTL, hierarchical TAD / PSYCHIC, TAD/ENCODE domains, and promoter annotations). It does **not** reimplement the linking logic.

Scientific parameters (promoter size = 600 bp, strand-aware definition, negative-coordinate filtering, ID mapping style, overlap conventions) are kept identical to Caitlin Mills’ original process (PLOS ONE 2020) and the modernized Python scripts.

## Quick start (for biologists)

1. Ensure static biocurated files are present in the repository root:
   - `CREbedDBenhancers_10092018`
   - `pantherGeneList.txt`
   - `tissuetable_10092018.txt`

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   # bedtools is also required (conda install -c bioconda bedtools)
   ```

3. Run environment checks:
   ```bash
   ./run_preprocess.sh check
   # or: python -m preprocess.cli check
   ```

4. Download public data and perform light cleaning:
   ```bash
   ./run_preprocess.sh all
   ```
   This places:
   - GTEx v7 `*signif_variant_gene_pairs.txt.gz` into `eQTL/GTEx_Analysis_v7_eQTL/`
   - PSYCHIC `*enh_1e-4.bed` into `heirarchicalTAD/data/`
   - ENCODE domain BEDs into `tad/`
   - Cleaned (negative-coordinate-removed) PSYCHIC files alongside the originals

5. Generate annotation intermediates (after obtaining Ensembl BioMart exports):
   ```bash
   ./run_preprocess.sh annotations \
       --gene-table TSS_default.txt \
       --promoter-out promoters_600bp.bed \
       --use-panther \
       --exon-table exon_file_default.txt \
       --exon-out eQTL/exons_genes.txt \
       --tss-table TSS_default.txt \
       --tss-out tad/TSSgenesbed
   ```

6. Run the original linking pipeline as usual:
   ```bash
   ./steps.sh /path/to/output
   ```

## What is produced

| Intermediate | Format | Produced by | Consumed by |
|--------------|--------|-------------|-------------|
| 600 bp promoter BED | `chr start end gene` (strand-aware: + = [TSS-600, TSS], − = [TSS, TSS+600]) | `annotations` | ChIA-PET / TAD promoter steps |
| eQTL BED-like | `chr start start variant_id gene_pval_tissue` | `eqtl` helper (or the existing `eqtl.py`) | Existing eQTL linking |
| Cleaned PSYCHIC BEDs | Original columns, negative coordinates removed | `psychic` | `heirarchicalTAD.py` |
| Formatted domain/TAD BEDs | chr-prefixed, no headers, valid intervals | `tad` helper | TAD scripts |

## Data sources

- **GTEx v7 eQTL**: https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz
- **PSYCHIC hierarchical TADs**: https://www.cs.huji.ac.il/~tommy/PSYCHIC/ (`*enh_1e-4.bed`)
- **ENCODE domains**: ENCFF274VJU, ENCFF588KUZ (nested contact domains, hg19)
- **Ensembl BioMart (GRCh37)**: gene/TSS and exon exports (see `docs/README_exon_file_default.txt.md` for exact attribute selections). Automated REST download is not required; place the TSV exports and run the annotations command.

## Environment checks

The `check` command verifies:
- Python ≥ 3.7
- `bedtools` on PATH
- Required Python packages (pandas, requests, pybedtools)
- Presence of the three static biocurated files

Clear error messages are raised if anything is missing.

## Design notes

- New code lives entirely under `preprocess/` and does not modify the existing linking scripts.
- Biocurated files remain static inputs.
- The workflow is intentionally simple for infrequent updates.
- Fidelity to the original 2020 process is prioritized over algorithmic novelty.

For questions about the original linking logic, refer to the READMEs in `eQTL/`, `heirarchicalTAD/`, `tad/`, and `chia_pet/`, and to the PLOS ONE paper.
