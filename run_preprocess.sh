#!/usr/bin/env bash
# Thin orchestration wrapper for the PEREGRINE pre-processing pipeline.
# Designed for biologists: simple commands, clear errors.

set -euo pipefail

REPO_DIR="$(CDPATH= cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_DIR"

PYTHON="${PYTHON:-python3}"

usage() {
    cat <<EOF
Usage: $0 <command> [options]

Commands:
  check                 Run environment and static-file checks
  download [steps]      Download public data (gtex,psychic,encode or subset)
  annotations ...       Build promoter/exon/TSS files (see --help)
  psychic               Clean PSYCHIC BEDs (remove negative coordinates)
  all                   download + full preprocessing (annotations, eQTL, TAD, PSYCHIC)

Examples:
  $0 check
  $0 download gtex,psychic
  $0 all
  $0 annotations --gene-table genes.tsv --promoter-out promoters.bed

Requires: Python 3.7+, bedtools, pandas, requests, pybedtools
Static inputs (must exist in repo root):
  CREbedDBenhancers_10092018
  pantherGeneList.txt
  tissuetable_10092018.txt
EOF
}

if [[ $# -lt 1 ]]; then
    usage
    exit 1
fi

CMD="$1"
shift || true

case "$CMD" in
    check)
        "$PYTHON" -m preprocess.cli check
        ;;
    download)
        "$PYTHON" -m preprocess.cli download --steps "${1:-gtex,psychic,encode}"
        ;;
    annotations|psychic|eqtl|tad)
        "$PYTHON" -m preprocess.cli "$CMD" "$@"
        ;;
    all)
        "$PYTHON" -m preprocess.cli all "$@"
        ;;
    help|--help|-h)
        usage
        ;;
    *)
        echo "Unknown command: $CMD" >&2
        usage
        exit 1
        ;;
esac
