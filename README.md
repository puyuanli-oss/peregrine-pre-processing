# PEREGRINE preprocessing — isolated pipeline

This repository contains **only preprocessing code**, with no downstream PEREGRINE linking algorithms and no bundled biological datasets. It downloads or accepts raw inputs and writes preprocessing artifacts under `data/processed/`. All `data/` contents are ignored by Git.

## Installation and tests

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python -m unittest -v tests_preprocessing
```

## Run

```bash
./run_preprocess.sh all
# Offline, with local GRCh37 GTF and pre-populated data/raw/{gtex,psychic,encode}:
./run_preprocess.sh all --offline --gtf /path/to/genes.gtf.gz
```

The raw input locations are `data/raw/gtex/`, `data/raw/psychic/`, and `data/raw/encode/`. Generated artifacts are in `data/processed/` (including `annotations/`, `eqtl/`, and `tad/`). The `all` command does not run downstream linking algorithms. `--steps` selects among `gtex,psychic,encode`; annotation preparation runs in all cases.

Optional PANTHER filtering (`--use-panther`) requires placing `pantherGeneList.txt` at repository root; it is **not bundled**. The curated enhancer resource `CREbedDBenhancers_10092018` is also **not bundled** and is not necessary for the preprocessing transformations in `all`; downstream linking may require it.

## Limitations

Offline tests validate synthetic transformations, not a complete production run. Network sources and genome assembly compatibility must be verified on real inputs. PSYCHIC download currently permits partial retrieval; review its warnings before treating outputs as complete. The GTF-derived TSS representation is gene-level and may differ from transcript-level legacy exports.
