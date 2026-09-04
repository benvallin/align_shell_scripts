#!/usr/bin/env bash
#
# merge_pb_fastq.sh
#
# Concatenates per-cell FASTQ files into pseudobulk FASTQ files, in parallel
# across pseudobulk groups, using a list of pseudobulk IDs and per-group
# per-read FASTQ path lists.
#
# Expected input file layout under LISTS_DIR:
#   pb_id_list_<LOG_ID>
#   <pb_id>_R1_fastq_list_<LOG_ID>
#   <pb_id>_R2_fastq_list_<LOG_ID>
#
# Output (written under an OUT_DIR/LOG_ID subdirectory, created if needed):
#   <OUT_DIR>/<LOG_ID>/<pb_id>_R1_fq.gz
#   <OUT_DIR>/<LOG_ID>/<pb_id>_R2_fq.gz
#
# Usage:
#   ./merge_pb_fastq.sh -l LISTS_DIR -i LOG_ID -o OUT_DIR [-P N] [-h]
#
# Options:
#   -l LISTS_DIR   Directory containing the pb_id_list and per-group fastq lists (required)
#   -i LOG_ID      Log/run identifier suffix used in the list filenames (required)
#   -o OUT_DIR     Output directory for merged pseudobulk FASTQ files (required)
#   -P N           Number of pseudobulk groups to process in parallel (default: 4)
#   -h             Show this help message and exit

set -euo pipefail

usage() {
  # Print only the leading header comment block (from line 2 up to the
  # first blank line), not every '#'-prefixed line in the script.
  sed -n '2,/^$/p' "$0" | sed -e '/^$/d' -e 's/^# \{0,1\}//'
  exit "${1:-0}"
}

# ---- Defaults ----

PARALLEL=4

# ---- Parse arguments ----

while getopts ":l:i:o:P:h" opt; do
  case "${opt}" in
    l) LISTS_DIR="${OPTARG}" ;;
    i) LOG_ID="${OPTARG}" ;;
    o) OUT_DIR="${OPTARG}" ;;
    P) PARALLEL="${OPTARG}" ;;
    h) usage 0 ;;
    \?) echo "Error: invalid option -${OPTARG}" >&2; usage 1 ;;
    :) echo "Error: option -${OPTARG} requires an argument" >&2; usage 1 ;;
  esac
done

# ---- Validate required arguments ----

: "${LISTS_DIR:?Error: -l LISTS_DIR is required}"
: "${LOG_ID:?Error: -i LOG_ID is required}"
: "${OUT_DIR:?Error: -o OUT_DIR is required}"

if [[ ! -d "${LISTS_DIR}" ]]; then
  echo "Error: LISTS_DIR does not exist: ${LISTS_DIR}" >&2
  exit 1
fi

PB_ID_LIST="${LISTS_DIR}/pb_id_list_${LOG_ID}"
if [[ ! -f "${PB_ID_LIST}" ]]; then
  echo "Error: pb_id_list file not found: ${PB_ID_LIST}" >&2
  exit 1
fi

if ! [[ "${PARALLEL}" =~ ^[0-9]+$ ]] || [[ "${PARALLEL}" -lt 1 ]]; then
  echo "Error: -P must be a positive integer (got '${PARALLEL}')" >&2
  exit 1
fi

FINAL_OUT_DIR="${OUT_DIR}/${LOG_ID}"
mkdir -p "${FINAL_OUT_DIR}"

echo "LISTS_DIR     : ${LISTS_DIR}"
echo "LOG_ID        : ${LOG_ID}"
echo "OUT_DIR       : ${OUT_DIR}"
echo "FINAL_OUT_DIR : ${FINAL_OUT_DIR}"
echo "PARALLEL      : ${PARALLEL}"
echo

# ---- Concatenation function ----

concat_pb() {
  local i="$1"
  local lists_dir="$2"
  local log_id="$3"
  local out_dir="$4"

  local r1_list="${lists_dir}/${i}_R1_fastq_list_${log_id}"
  local r2_list="${lists_dir}/${i}_R2_fastq_list_${log_id}"

  if [[ ! -f "${r1_list}" ]]; then
    echo "Error [${i}]: missing R1 list: ${r1_list}" >&2
    return 1
  fi
  if [[ ! -f "${r2_list}" ]]; then
    echo "Error [${i}]: missing R2 list: ${r2_list}" >&2
    return 1
  fi

  echo "Processing pseudobulk ${i}..."
  xargs cat < "${r1_list}" > "${out_dir}/${i}_R1_fq.gz"
  xargs cat < "${r2_list}" > "${out_dir}/${i}_R2_fq.gz"
}
export -f concat_pb

# ---- Run concatenation in parallel across pseudobulk groups ----

xargs -P "${PARALLEL}" -I{} bash -c 'concat_pb "$@"' _ {} "${LISTS_DIR}" "${LOG_ID}" "${FINAL_OUT_DIR}" \
  < "${PB_ID_LIST}"

# ---- Integrity check ----

echo
echo "Running integrity checks..."
FAILED=0
while read -r i; do
  [[ -z "${i}" ]] && continue
  for r in R1 R2; do
    f="${FINAL_OUT_DIR}/${i}_${r}_fq.gz"
    if [[ -f "${f}" ]] && gzip -t "${f}" 2>/dev/null; then
      echo "OK: ${f}"
    else
      echo "CORRUPT or MISSING: ${f}" >&2
      FAILED=1
    fi
  done
done < "${PB_ID_LIST}"

if [[ "${FAILED}" -ne 0 ]]; then
  echo
  echo "One or more output files failed integrity checks. See above." >&2
  exit 1
fi

echo
echo "All pseudobulk FASTQ files merged and verified successfully."