"""Command-line interface for the PEREGRINE pre-processing pipeline.

Usage examples:
    python -m preprocess.cli check
    python -m preprocess.cli download --steps gtex,psychic,encode
    python -m preprocess.cli annotations --gene-table genes.tsv --promoter-out promoters.bed
    python -m preprocess.cli all
"""

import argparse
import sys
from pathlib import Path

from .checks import run_all_checks
from .download import download_gtex_v7, download_psychic, download_encode_domains
from .annotations import generate_promoter_bed, process_exons, process_tss_bed
from .psychic import process_psychic_directory, remove_negative_coordinates
from .eqtl import eqtl_to_bedlike
from .tad import format_domain_bed
from .auto import prepare_annotations, prepare_eqtl, prepare_tad


def repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def cmd_check(args):
    run_all_checks(repo_root())


def cmd_download(args):
    root = repo_root()
    steps = set(args.steps.split(",")) if args.steps else {"gtex", "psychic", "encode"}
    if "gtex" in steps:
        download_gtex_v7(root / "data" / "raw" / "gtex", force=args.force)
    if "psychic" in steps:
        download_psychic(root / "data" / "raw" / "psychic", force=args.force)
    if "encode" in steps:
        download_encode_domains(root / "data" / "raw" / "encode", force=args.force)


def cmd_annotations(args):
    root = repo_root()
    if args.gene_table and args.promoter_out:
        generate_promoter_bed(
            Path(args.gene_table),
            Path(args.promoter_out),
            promoter_size=args.promoter_size,
            panther_list=root / "pantherGeneList.txt" if args.use_panther else None,
        )
    if args.exon_table and args.exon_out:
        process_exons(Path(args.exon_table), Path(args.exon_out), args.mapping)
    if args.tss_table and args.tss_out:
        process_tss_bed(Path(args.tss_table), Path(args.tss_out), args.mapping)


def cmd_psychic(args):
    if args.input and args.output:
        remove_negative_coordinates(Path(args.input), Path(args.output), args.tossed)
    elif args.dir:
        process_psychic_directory(Path(args.dir), Path(args.outdir or args.dir))
    else:
        # Default: process the data directory if present
        data_dir = repo_root() / "heirarchicalTAD" / "data"
        if data_dir.is_dir():
            process_psychic_directory(data_dir, data_dir)
        else:
            print("No PSYCHIC data directory found. Download first or specify --dir.")


def cmd_eqtl(args):
    if not args.input or not args.output:
        print("Usage: eqtl --input VARIANT_GENE_FILE --output BEDLIKE")
        return
    eqtl_to_bedlike(Path(args.input), Path(args.output), args.tissue)


def cmd_tad(args):
    if not args.input or not args.output:
        print("Usage: tad --input RAW.bed --output CLEAN.bed")
        return
    format_domain_bed(Path(args.input), Path(args.output))


def cmd_all(args):
    """Run selected acquisition and every applicable transformation, failing on missing inputs."""
    root = repo_root()
    steps = set(args.steps.split(',')) if args.steps else {'gtex', 'psychic', 'encode'}
    invalid = steps - {'gtex', 'psychic', 'encode'}
    if invalid:
        raise ValueError(f'Unknown download steps: {sorted(invalid)}')
    if not args.offline:
        cmd_download(args)
    if 'psychic' in steps:
        data_dir = root / 'data' / 'raw' / 'psychic'
        beds = list(data_dir.glob('*.bed')) if data_dir.exists() else []
        if not beds:
            raise FileNotFoundError('No PSYCHIC BED files; download or supply them')
        process_psychic_directory(data_dir, data_dir)
    outputs = prepare_annotations(root, gtf=args.gtf, download=not args.offline, force=args.force,
                                  use_panther=args.use_panther)
    if 'gtex' in steps:
        print('eQTL files converted:', prepare_eqtl(root))
    if 'encode' in steps:
        print('TAD files formatted:', prepare_tad(root))
    print('Completed preprocessing:', {k: str(v) for k, v in outputs.items()})


def main():
    parser = argparse.ArgumentParser(
        description="PEREGRINE pre-processing pipeline (annotations, eQTL, PSYCHIC, domains)"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_check = sub.add_parser("check", help="Run environment and static-file checks")
    p_check.set_defaults(func=cmd_check)

    p_dl = sub.add_parser("download", help="Download public data")
    p_dl.add_argument("--steps", default="gtex,psychic,encode",
                      help="Comma-separated: gtex,psychic,encode")
    p_dl.add_argument("--force", action="store_true")
    p_dl.set_defaults(func=cmd_download)

    p_ann = sub.add_parser("annotations", help="Build promoter / exon / TSS intermediates")
    p_ann.add_argument("--gene-table", help="Ensembl BioMart gene/TSS TSV")
    p_ann.add_argument("--promoter-out", help="Output 600 bp promoter BED")
    p_ann.add_argument("--promoter-size", type=int, default=600)
    p_ann.add_argument("--use-panther", action="store_true",
                       help="Filter promoters to genes present in pantherGeneList.txt")
    p_ann.add_argument("--exon-table", help="Ensembl BioMart exon TSV")
    p_ann.add_argument("--exon-out", help="Output exons_genes.txt")
    p_ann.add_argument("--tss-table", help="Ensembl BioMart TSS TSV")
    p_ann.add_argument("--tss-out", help="Output TSSgenesbed")
    p_ann.add_argument("--mapping", help="Optional patch chromosome mapping file")
    p_ann.set_defaults(func=cmd_annotations)

    p_psy = sub.add_parser("psychic", help="Remove negative coordinates from PSYCHIC BEDs")
    p_psy.add_argument("--input")
    p_psy.add_argument("--output")
    p_psy.add_argument("--tossed")
    p_psy.add_argument("--dir")
    p_psy.add_argument("--outdir")
    p_psy.set_defaults(func=cmd_psychic)

    p_eq = sub.add_parser("eqtl", help="Convert GTEx pairs file to BED-like intermediate")
    p_eq.add_argument("--input")
    p_eq.add_argument("--output")
    p_eq.add_argument("--tissue")
    p_eq.set_defaults(func=cmd_eqtl)

    p_tad = sub.add_parser("tad", help="Format ENCODE domain/TAD BED")
    p_tad.add_argument("--input")
    p_tad.add_argument("--output")
    p_tad.set_defaults(func=cmd_tad)

    p_all = sub.add_parser("all", help="Download and preprocess annotations, GTEx, PSYCHIC, and ENCODE")
    p_all.add_argument("--steps", default="gtex,psychic,encode")
    p_all.add_argument("--force", action="store_true")
    p_all.add_argument("--offline", action="store_true", help="Use locally supplied data without network")
    p_all.add_argument("--gtf", help="Local GRCh37 Ensembl GTF or GTF.gz (otherwise download release 87)")
    p_all.add_argument("--use-panther", action="store_true")
    p_all.set_defaults(func=cmd_all)

    args = parser.parse_args()
    try:
        args.func(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
