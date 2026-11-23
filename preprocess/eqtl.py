"""eQTL pre-processing helpers.

The existing eQTL/eqtl.py expects:
- *signif_variant_gene_pairs.txt.gz in eQTL/GTEx_Analysis_v7_eQTL/
- exons_genes.txt in eQTL/

This module prepares the BED-like intermediate that eqtl_process produces,
so that biologists can inspect or reuse it independently.
"""

from pathlib import Path
import gzip
import pandas as pd


def eqtl_to_bedlike(variant_gene_file: Path, output_bed: Path, tissue: str = None):
    """Convert a GTEx signif_variant_gene_pairs file to the BED-like format
    produced by eQTL/eqtl.py::eqtl_process.

    Output columns (tab-separated, no header):
    chrom  start  start  variant_id  gene_id_pval_tissue
    """
    variant_gene_file = Path(variant_gene_file)
    output_bed = Path(output_bed)
    output_bed.parent.mkdir(parents=True, exist_ok=True)

    opener = gzip.open if variant_gene_file.suffix == ".gz" else open
    mode = "rt" if variant_gene_file.suffix == ".gz" else "r"

    if tissue is None:
        tissue = variant_gene_file.name.split(".")[0]

    records = []
    with opener(variant_gene_file, mode) as f:
        header = f.readline().rstrip("\n").split("\t")
        # Expected: variant_id, gene_id, pval_nominal, ...
        try:
            vid_idx = header.index("variant_id")
            gid_idx = header.index("gene_id")
            pval_idx = header.index("pval_nominal")
        except ValueError as exc:
            raise ValueError(f"{variant_gene_file}: missing required columns") from exc

        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(vid_idx, gid_idx, pval_idx):
                continue
            variant_id = parts[vid_idx]
            gene_id = parts[gid_idx]
            pval = parts[pval_idx]
            # variant_id format: 1_100169983_A_G_b37
            vparts = variant_id.split("_")
            if len(vparts) < 2:
                continue
            chrom = "chr" + vparts[0]
            start = vparts[1]
            key = f"{chrom}\t{start}\t{start}\t{variant_id}\t{gene_id}_{pval}_{tissue}"
            records.append(key)

    with open(output_bed, "w") as out:
        out.write("\n".join(records) + ("\n" if records else ""))
    print(f"Wrote {len(records)} eQTL intervals to {output_bed}")
    return output_bed
