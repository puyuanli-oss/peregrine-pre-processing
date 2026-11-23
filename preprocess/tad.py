"""ENCODE domain / TAD file formatting helpers."""

from pathlib import Path


def format_domain_bed(input_bed: Path, output_bed: Path):
    """Ensure domain/TAD BED is clean: chr-prefixed, no header, tab-separated,
    and free of negative coordinates.
    """
    input_bed = Path(input_bed)
    output_bed = Path(output_bed)
    output_bed.parent.mkdir(parents=True, exist_ok=True)

    kept = 0
    with open(input_bed) as src, open(output_bed, "w") as out:
        for line in src:
            if not line.strip() or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            if start < 0 or end < 0 or end <= start:
                continue
            rest = parts[3:] if len(parts) > 3 else []
            out.write("\t".join([chrom, str(start), str(end)] + rest) + "\n")
            kept += 1
    print(f"Formatted {kept} domain intervals -> {output_bed}")
    return output_bed
