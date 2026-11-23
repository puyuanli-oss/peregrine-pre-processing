"""Streaming GRCh37 GTF annotation acquisition and preprocessing orchestration."""
import gzip
import re
from pathlib import Path
from urllib.request import urlopen
from .annotations import generate_promoter_bed, process_exons, process_tss_bed
from .eqtl import eqtl_to_bedlike
from .tad import format_domain_bed

# Explicitly pinned GRCh37 release; override URL for mirrors or local cached input.
GTF_URL = 'https://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz'


def fetch_gtf(destination, url=GTF_URL, force=False):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size and not force:
        return destination
    temp = destination.with_name(destination.name + '.partial')
    try:
        with urlopen(url, timeout=90) as src, open(temp, 'wb') as dst:
            while True:
                chunk = src.read(1024 * 1024)
                if not chunk:
                    break
                dst.write(chunk)
        with gzip.open(temp, 'rb') as f:
            if not f.read(1):
                raise ValueError('Downloaded GTF is empty')
        temp.replace(destination)
    finally:
        temp.unlink(missing_ok=True)
    return destination


def gtf_to_tables(gtf, directory):
    """Generate BioMart-compatible tables using 1-based Ensembl GTF coordinates.

    Promoter/TSS use gene feature start/end; exon records are retained per GTF row.
    Genome assembly must be GRCh37; this function cannot infer assembly from arbitrary files.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    gene_table = directory / 'genes_from_gtf.tsv'
    exon_table = directory / 'exons_from_gtf.tsv'
    opener = gzip.open if str(gtf).endswith('.gz') else open
    gene_count = exon_count = 0
    with opener(gtf, 'rt') as source, open(gene_table, 'w') as genes, open(exon_table, 'w') as exons:
        genes.write('Gene stable ID\tChromosome/scaffold name\tTranscription start site (TSS)\tStrand\n')
        exons.write('Gene stable ID\tChromosome/scaffold name\tExon region start (bp)\tExon region end (bp)\tGene type\n')
        for line in source:
            if not line or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9 or fields[2] not in ('gene', 'exon'):
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = fields
            if not re.fullmatch(r'(?:[0-9]{1,2}|X|Y|MT)', chrom):
                continue
            match = re.search(r'(?:^|;\s*)gene_id "([^"]+)"', attrs)
            if not match or strand not in ('+', '-'):
                continue
            gene = match.group(1).split('.')[0]
            if feature == 'gene':
                tss = start if strand == '+' else end
                genes.write(f'{gene}\t{chrom}\t{tss}\t{strand}\n')
                gene_count += 1
            else:
                kind = re.search(r'(?:^|;\s*)gene_biotype "([^"]+)"', attrs)
                kind = kind or re.search(r'(?:^|;\s*)gene_type "([^"]+)"', attrs)
                exons.write(f'{gene}\t{chrom}\t{start}\t{end}\t{kind.group(1) if kind else "unknown"}\n')
                exon_count += 1
    if not gene_count or not exon_count:
        raise ValueError(f'GTF lacked gene or exon features: genes={gene_count}, exons={exon_count}')
    return gene_table, exon_table


def prepare_annotations(root, gtf=None, download=True, force=False, use_panther=False):
    root = Path(root)
    inputs = root / 'data' / 'processed' / 'annotations'
    if gtf is None:
        gtf = root / 'data' / 'processed' / 'sources' / 'Homo_sapiens.GRCh37.87.gtf.gz'
        if download:
            fetch_gtf(gtf, force=force)
        elif not gtf.exists():
            raise FileNotFoundError(f'Provide --gtf or download the GRCh37 GTF: {gtf}')
    genes, exons = gtf_to_tables(gtf, inputs)
    outputs = {
        'promoters': root / 'data' / 'processed' / 'promoters_600bp.bed',
        'exons': root / 'data' / 'processed' / 'annotations' / 'exons_genes.txt',
        'tss': root / 'data' / 'processed' / 'annotations' / 'TSSgenesbed',
    }
    generate_promoter_bed(genes, outputs['promoters'], panther_list=root / 'pantherGeneList.txt' if use_panther else None)
    process_exons(exons, outputs['exons'])
    process_tss_bed(genes, outputs['tss'])
    for name, path in outputs.items():
        if not path.exists() or path.stat().st_size == 0:
            raise RuntimeError(f'Empty {name} output: {path}')
    return outputs


def prepare_eqtl(root):
    root = Path(root)
    files = sorted((root / 'data' / 'raw' / 'gtex' / 'GTEx_Analysis_v7_eQTL').glob('*signif_variant_gene_pairs.txt.gz'))
    if not files:
        raise FileNotFoundError('No GTEx significant variant-gene pairs found')
    outdir = root / 'data' / 'processed' / 'eqtl'
    for file in files:
        eqtl_to_bedlike(file, outdir / (file.name[:-7] + '.bed'), tissue=file.name.split('.signif_')[0])
    return len(files)


def prepare_tad(root):
    root = Path(root)
    files = sorted((root / 'data' / 'raw' / 'encode').glob('ENCFF*.bed'))
    if not files:
        raise FileNotFoundError('No ENCODE domain BED files found')
    outdir = root / 'data' / 'processed' / 'tad'
    for file in files:
        format_domain_bed(file, outdir / file.name)
    return len(files)
