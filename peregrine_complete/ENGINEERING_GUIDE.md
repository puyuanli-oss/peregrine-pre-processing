# PEREGRINE engineering foundation

## Status and scope

This is an **engineering foundation**, not a claim of production biological validation. It preserves the existing linking algorithms and extends the preprocessing package with a supervised orchestration entry point, audit logs, ChIA-PET preflight, and final-link structural checks. ChIA-PET linking already exists in `chia_pet/ChIA_PET.py`; this release integrates it into the single top-level run. The input format is **post-Mango** eight-column interaction data, not raw FASTQ or BAM; mapping, peak calling and Mango itself are **not** automated.

## Installation and commands

Install the requirements in `requirements.txt`, Python 3.9+, and BEDTools. The legacy full linker requires a Python virtual environment at `venv/bin/python`; consult the original READMEs for its dependencies. Use the original curated `CREbedDBenhancers_10092018` in the repository root. It is intentionally not fabricated or bundled. Provide real Mango sample files at paths specified by `chia_pet/samples.tsv`, or supply `--manifest` with a different TSV. The default gene table is `chia_pet/hg19genes.txt`.

```bash
python -m unittest -v tests_preprocessing tests_engineering
python peregrine_engineering.py preflight --output results
python peregrine_engineering.py preprocess --output results
python peregrine_engineering.py run --output results
python peregrine_engineering.py validate --output results
```

Use `--offline --gtf /path/to/GRCh37.gtf.gz` to avoid annotation downloads when the other downloaded inputs are already cached. `--skip-preprocess` runs linking on existing preprocessed inputs. `--manifest` and `--genes` select the ChIA sources. **No command automatically downloads Mango ChIA data**; their availability, license and provenance must be established by the operator.

## Output contract

After a successful full run, `results/engineering/report.json` records command status, stages, log paths, SHA-256 hashes and row counts for the following files:

- `results/eQTL/linksDBnumeqtl`
- `results/heirarchicalTAD/PSYCHIClinksDB`
- `results/chia_pet/linksDBchia`
- `results/tad/linksDBtad`

Missing outputs fail validation. Empty existing outputs are permitted because a valid negative result can have zero links. The validator checks four nonempty tab-separated fields for each nonempty row, and ChIA assay `1`. It does **not** prove statistical validity, reference-build agreement, completeness of inputs, biological equivalence, or correctness of every link. Each run creates stage logs; reruns overwrite the same log names and report. For immutable provenance, archive the engineering folder per run.

## Critical limitations and next validation gates

1. Obtain the curated enhancer file and post-Mango ChIA datasets from authoritative sources, with explicit assembly and version metadata.
2. Confirm Ensembl GTF-derived gene-level annotations are biologically equivalent to the original BioMart exports. Transcript-level TSS differences remain a known risk.
3. Validate chromosome naming, coordinate conventions, GRCh37/hg19 assembly and tissue-code mapping across all inputs.
4. Run on a production-sized dataset with enough disk/RAM/CPU and a pinned environment; compare output link sets against a trusted PEREGRINE baseline.
5. Inspect complete stage logs and audit output; verify scientific output distributions and test representative enhancer-gene links manually.
6. Add integration tests against a distributable, licensed fixture covering all four legacy linking algorithms. Current offline tests cover transformations and engineering contracts, not full real-data execution.

**Safety:** Never interpret a passing structural validator as a scientific quality certificate.
