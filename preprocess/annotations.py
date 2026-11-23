"""Annotation processing: 600 bp strand-aware promoters and exon files.

Scientific parameters match Caitlin Mills' original PEREGRINE process
and the modernized ChIA_PET.py / TSSgenesbed_generation.py / exon_processing.py.
"""

from pathlib import Path
import pandas as pd


def load_patch_mapping(mapping_file: Path) -> dict:
    """Load patch chromosome mapping (from biostars/Ensembl style)."""
    mapping = {}
    if not mapping_file or not Path(mapping_file).is_file():
        return mapping
    with open(mapping_file) as f:
        for line in f:
            parts = [p for p in line.strip().split() if p]
            if len(parts) >= 2:
                mapping[parts[0]] = "chr" + parts[1]
    return mapping


def normalize_chromosome(chrom: str, mapping: dict) -> str:
    chrom = str(chrom).strip()
    if chrom in mapping:
        return mapping[chrom]
    if chrom.startswith("chr"):
        return chrom
    return "chr" + chrom


def generate_promoter_bed(
    gene_table: Path,
    output_bed: Path,
    promoter_size: int = 600,
    panther_list: Path = None,
):
    """Generate strand-aware 600 bp promoter BED.

    Matches ChIA_PET.py promoter_lines logic:
    - + strand (1): [TSS - 600, TSS]
    - - strand (-1): [TSS, TSS + 600]
    - Chromosome without leading 'chr' in input; output with 'chr'
    - Last record per gene wins if duplicates; only genes in PANTHER if list given
    """
    df = pd.read_csv(gene_table, sep="\t", dtype=str)
    # Expected columns (flexible matching)
    col_map = {
        "gene": ["Gene stable ID", "gene_id", "ensembl_gene_id"],
        "tss": ["Transcription start site (TSS)", "tss", "TSS"],
        "strand": ["Strand", "strand"],
        "chrom": ["Chromosome/scaffold name", "chromosome", "chrom"],
    }
    resolved = {}
    for key, candidates in col_map.items():
        for c in candidates:
            if c in df.columns:
                resolved[key] = c
                break
        else:
            raise ValueError(
                f"Could not find column for '{key}' in {gene_table}. "
                f"Columns present: {list(df.columns)}"
            )

    panther_genes = None
    if panther_list and Path(panther_list).is_file():
        panther_genes = set()
        with open(panther_list) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    panther_genes.add(parts[1].split(".")[0])

    genes = {}
    for _, row in df.iterrows():
        gene = str(row[resolved["gene"]]).split(".")[0]
        if panther_genes is not None and gene not in panther_genes:
            continue
        strand = str(row[resolved["strand"]]).strip()
        if strand not in ("1", "-1", "+", "-"):
            continue
        try:
            tss = int(float(row[resolved["tss"]]))
        except (ValueError, TypeError):
            continue
        chrom = str(row[resolved["chrom"]]).strip()
        if chrom.startswith("chr"):
            chrom_out = chrom
        else:
            chrom_out = "chr" + chrom

        if strand in ("1", "+"):
            start, end = tss - promoter_size, tss
        else:
            start, end = tss, tss + promoter_size

        if start < 0:
            start = 0
        # Last eligible record wins (as in original)
        genes[gene] = f"{chrom_out}\t{start}\t{end}\t{gene}"

    output_bed = Path(output_bed)
    output_bed.parent.mkdir(parents=True, exist_ok=True)
    records = sorted(set(genes.values()))
    with open(output_bed, "w") as out:
        for rec in records:
            out.write(rec + "\n")
    print(f"Wrote {len(records)} promoter intervals to {output_bed}")
    return output_bed


def process_exons(
    exon_table: Path,
    output_file: Path,
    mapping_file: Path = None,
):
    """Convert Ensembl BioMart exon export to the exons_genes.txt format
    expected by eQTL/eqtl.py eliminateOverlap.
    """
    df = pd.read_csv(exon_table, sep="\t", dtype=str)
    mapping = load_patch_mapping(mapping_file) if mapping_file else {}

    col_map = {
        "chrom": ["Chromosome/scaffold name", "chromosome"],
        "start": ["Exon region start (bp)", "exon_start"],
        "end": ["Exon region end (bp)", "exon_end"],
        "gene": ["Gene stable ID", "gene_id"],
        "gene_type": ["Gene type", "gene_biotype"],
    }
    resolved = {}
    for key, candidates in col_map.items():
        for c in candidates:
            if c in df.columns:
                resolved[key] = c
                break
        else:
            if key == "gene_type":
                continue
            raise ValueError(
                f"Missing column for '{key}' in {exon_table}. Columns: {list(df.columns)}"
            )

    df = df.copy()
    df["chrom_out"] = df[resolved["chrom"]].apply(lambda x: normalize_chromosome(x, mapping))
    df["start_out"] = df[resolved["start"]]
    df["end_out"] = df[resolved["end"]]
    df["blank"] = "-"
    df["gene_out"] = df[resolved["gene"]].str.split(".").str[0]
    df["type_out"] = df[resolved.get("gene_type", resolved["gene"])] if "gene_type" in resolved else "protein_coding"
    df["known"] = "KNOWN"
    df["lciz"] = "LCIZ"

    out_cols = ["chrom_out", "start_out", "end_out", "blank", "gene_out", "type_out", "known", "lciz"]
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    df[out_cols].to_csv(output_file, sep="\t", header=False, index=False)
    print(f"Wrote {len(df)} exon records to {output_file}")
    return output_file


def process_tss_bed(
    tss_table: Path,
    output_file: Path,
    mapping_file: Path = None,
):
    """Convert Ensembl BioMart TSS export to the simple TSS BED used by some TAD steps.
    Format: chrom\tstart\tend\tgene  (point intervals at TSS)
    """
    df = pd.read_csv(tss_table, sep="\t", dtype=str)
    mapping = load_patch_mapping(mapping_file) if mapping_file else {}

    chrom_col = next((c for c in ["Chromosome/scaffold name", "chromosome"] if c in df.columns), None)
    tss_col = next((c for c in ["Transcription start site (TSS)", "tss"] if c in df.columns), None)
    gene_col = next((c for c in ["Gene stable ID", "gene_id"] if c in df.columns), None)
    if not all([chrom_col, tss_col, gene_col]):
        raise ValueError(f"Could not identify required columns in {tss_table}")

    df = df.copy()
    df["chrom_out"] = df[chrom_col].apply(lambda x: normalize_chromosome(x, mapping))
    df["tss"] = df[tss_col]
    df["gene"] = df[gene_col].str.split(".").str[0]

    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    df[["chrom_out", "tss", "tss", "gene"]].to_csv(output_file, sep="\t", header=False, index=False)
    print(f"Wrote {len(df)} TSS records to {output_file}")
    return output_file
