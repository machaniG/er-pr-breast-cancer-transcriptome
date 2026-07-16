#!/usr/bin/env bash
set -euo pipefail

# ==============================================================
# qc_raw_reads.sh
# Raw Read Quality Control with FastQC + MultiQC
# ==============================================================

# Usage:
#   bash qc_raw_reads.sh
#
# Optional:
#   edit FASTQ_DIR / RESULTS_DIR below if your files live elsewhere.
#
# Expected:
#   FASTQ files with extension .fastq.gz
# ==============================================================

# =========================
# CONFIGURATION
# =========================
FASTQ_DIR="./raw-data"
RESULTS_DIR="./qc-output"
FASTQC_DIR="${RESULTS_DIR}/fastqc"
MULTIQC_DIR="${RESULTS_DIR}/multiqc"

# macOS: sysctl, Linux: nproc
CPU_COUNT="$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)"
THREADS=$(( CPU_COUNT > 8 ? 8 : CPU_COUNT ))

# =========================
# COLORS
# =========================
BOLD='\033[1m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
RESET='\033[0m'

# =========================
# HELPERS
# =========================
section() {
  printf "\n${BOLD}── %s ──────────────────────────────────────────${RESET}\n" "$1"
}

ok() {
  printf " ${GREEN}[OK]${RESET} %s\n" "$*"
}

info() {
  printf " ${YELLOW}[!!]${RESET} %s\n" "$*"
}

die() {
  printf "\n ${RED}[ERR]${RESET} %s\n\n" "$*" >&2
  exit 1
}

check_tool() {
  local tool="$1"
  local install_hint="$2"

  if ! command -v "$tool" >/dev/null 2>&1; then
    die "'$tool' is not installed or not in PATH.
Install hint: $install_hint
Conda hint: conda install -c bioconda $tool"
  fi

  local tool_path tool_version
  tool_path="$(command -v "$tool")"
  tool_version="$("$tool" --version 2>&1 | head -n 1 || true)"

  ok "$tool"
  printf " path    : %s\n" "$tool_path"
  printf " version : %s\n" "$tool_version"
}

find_fastq_files() {
  [[ -d "$FASTQ_DIR" ]] || die "FASTQ folder does not exist: $FASTQ_DIR"

  shopt -s nullglob
  FASTQ_FILES=( "$FASTQ_DIR"/*.fastq.gz )
  FASTA_FILES=( "$FASTQ_DIR"/*.fasta.gz )
  shopt -u nullglob

  local count="${#FASTQ_FILES[@]}"
  [[ "$count" -gt 0 ]] || die "No .fastq.gz files found in: $FASTQ_DIR"

  printf "\n Found %d .fastq.gz file(s) in:\n" "$count"
  printf " %s\n" "$FASTQ_DIR"

  local f size
  for f in "${FASTQ_FILES[@]}"; do
    size="$(du -sh "$f" 2>/dev/null | cut -f1)"
    printf " %-50s %s\n" "$(basename "$f")" "[$size]"
  done

  if [[ "${#FASTA_FILES[@]}" -gt 0 ]]; then
    info "Found ${#FASTA_FILES[@]} .fasta.gz file(s) — skipping them."
    info "FASTA files have no quality scores, so FastQC cannot process them."
  fi
}

run_fastqc() {
  mkdir -p "$FASTQC_DIR"

  printf "\n Files   : %d total\n" "${#FASTQ_FILES[@]}"
  printf " Threads : %d\n" "$THREADS"
  printf " Output  : %s\n\n" "$FASTQC_DIR"

  fastqc "${FASTQ_FILES[@]}" \
    --outdir "$FASTQC_DIR" \
    --threads "$THREADS" \
    --extract

  ok "FastQC complete."
}

run_multiqc() {
  mkdir -p "$MULTIQC_DIR"

  printf "\n Scanning : %s\n" "$FASTQC_DIR"
  printf " Report   : %s/multiqc_report.html\n\n" "$MULTIQC_DIR"

  multiqc "$FASTQC_DIR" \
    --outdir "$MULTIQC_DIR" \
    --filename "multiqc_report" \
    --force \
    --quiet

  ok "MultiQC complete."
}

show_summary() {
  printf "\n${BOLD}══════════════════════════════════════════════════════════${RESET}\n"
  printf "${BOLD} All done!${RESET}\n"
  printf "${BOLD}══════════════════════════════════════════════════════════${RESET}\n\n"

  printf " Individual reports → %s/\n" "$FASTQC_DIR"
  printf " Summary report     → ${BOLD}%s/multiqc_report.html${RESET}\n\n" "$MULTIQC_DIR"

  printf "${BOLD} What to look for in MultiQC:${RESET}\n\n"
  printf " ${GREEN}[OK]${RESET} Per-base sequence quality should mostly be high.\n"
  printf " ${GREEN}[OK]${RESET} Adapter content should be near 0%%.\n"
  printf " ${GREEN}[OK]${RESET} GC content should look plausible for your organism.\n"
  printf " ${YELLOW}[!!]${RESET} Per-base sequence content warnings are common in RNA-seq.\n"
  printf " ${YELLOW}[!!]${RESET} Duplication warnings are also common in RNA-seq.\n\n"

  printf " Open the report manually:\n"
  printf " ${BOLD}open %s/multiqc_report.html${RESET}\n\n" "$MULTIQC_DIR"
  printf "${BOLD}══════════════════════════════════════════════════════════${RESET}\n\n"
}

# =========================
# MAIN
# =========================
printf "\n${BOLD}══════════════════════════════════════════════════════════${RESET}\n"
printf "${BOLD} Step 2: Raw Read Quality Control${RESET}\n"
printf " Directory : %s\n" "$(pwd)"
printf " Started   : %s\n" "$(date '+%Y-%m-%d %H:%M:%S')"
printf "${BOLD}══════════════════════════════════════════════════════════${RESET}\n"

section "1/4 Checking required tools"
check_tool "fastqc" "brew install fastqc"
check_tool "multiqc" "pip install multiqc"

section "2/4 Locating FASTQ files"
find_fastq_files

section "3/4 Running FastQC"
run_fastqc

section "4/4 Running MultiQC"
run_multiqc

show_summary

printf " Open the MultiQC report in your browser now? [Y/n]: "
read -r user_answer < /dev/tty 2>/dev/null || user_answer="y"

case "${user_answer:-y}" in
  [Nn]*)
    printf " Skipped. Open it anytime with:\n"
    printf " ${BOLD}open %s/multiqc_report.html${RESET}\n\n" "$MULTIQC_DIR"
    ;;
  *)
    printf " Opening...\n\n"
    open "$MULTIQC_DIR/multiqc_report.html"
    ;;
esac