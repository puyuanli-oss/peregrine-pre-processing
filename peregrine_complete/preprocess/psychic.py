"""PSYCHIC hierarchical TAD pre-processing: negative coordinate removal."""

from pathlib import Path


def remove_negative_coordinates(input_bed: Path, output_bed: Path, tossed_bed: Path = None):
    """Remove intervals with negative start or end coordinates.

    Matches the original nonegvalues.pl / heirarchicalTAD.py nonegvalues logic.
    """
    input_bed = Path(input_bed)
    output_bed = Path(output_bed)
    output_bed.parent.mkdir(parents=True, exist_ok=True)
    if tossed_bed:
        tossed_bed = Path(tossed_bed)
        tossed_bed.parent.mkdir(parents=True, exist_ok=True)

    kept = 0
    tossed = 0
    with open(input_bed) as src, open(output_bed, "w") as out:
        tossed_handle = open(tossed_bed, "w") if tossed_bed else None
        try:
            for line in src:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    if tossed_handle:
                        tossed_handle.write(line)
                    tossed += 1
                    continue
                if start < 0 or end < 0:
                    if tossed_handle:
                        tossed_handle.write(line)
                    tossed += 1
                else:
                    out.write(line if line.endswith("\n") else line + "\n")
                    kept += 1
        finally:
            if tossed_handle:
                tossed_handle.close()

    print(f"{input_bed.name}: kept {kept}, tossed {tossed} (negative coords)")
    return output_bed


def process_psychic_directory(input_dir: Path, output_dir: Path):
    """Process all *enh_1e-4.bed files, writing cleaned versions."""
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    beds = list(input_dir.glob("*.bed"))
    if not beds:
        raise FileNotFoundError(f"No BED files found in {input_dir}")
    results = []
    for bed in beds:
        out = output_dir / ("nn_" + bed.name)
        tossed = output_dir / ("tossed_" + bed.name)
        remove_negative_coordinates(bed, out, tossed)
        results.append(out)
    print(f"Processed {len(results)} PSYCHIC files into {output_dir}")
    return results
