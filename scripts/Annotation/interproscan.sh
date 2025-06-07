#!/usr/bin/env bash
#SBATCH -J ips-argania
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=160:00:00
#SBATCH --constraint=cal
#SBATCH --output=logs/interproscan.%J.out
#SBATCH --error=logs/interproscan.%J.err

# -------------------------------
# Load Environment
# -------------------------------
module purge
module load interproscan/5.67-99.0

set -e
trap 'echo "InterProScan failed at line $LINENO"; exit 1' ERR

# -------------------------------
# Define Paths
# -------------------------------
BASE="/mnt/home/users/pab_001_uma/nhabib/Argania-spinosa-sp"
INPUT="$BASE/transdecoder-output/transcripts_tpm_filtered.cleaned.pep"
OUTDIR="$BASE/interproscan-results"
IPRSCAN_TEMP="$BASE/temp"
mkdir -p "$OUTDIR" "$IPRSCAN_TEMP" logs

# -------------------------------
# Output Files with Timestamp
# -------------------------------
TS=$(date +%Y%m%d_%H%M%S)
OUTPUT_TSV="$OUTDIR/ips_output_${TS}.tsv"
OUTPUT_GFF="$OUTDIR/ips_output_${TS}.gff3"
LOGFILE="$OUTDIR/ips_log_${TS}.txt"

# -------------------------------
# Run InterProScan
# -------------------------------
echo "Running InterProScan on: $INPUT"
echo "Start: $(date)"
echo "Output: $OUTPUT_TSV"

interproscan.sh \
  -i "$INPUT" \
  -f tsv,gff3 \
  --goterms \
  --pathways \
  --iprlookup \
  --disable-precalc \
  --cpu "$SLURM_CPUS_PER_TASK" \
  -T "$IPRSCAN_TEMP" \
  -o "$OUTPUT_TSV"

echo "InterProScan completed"
echo "TSV Output: $OUTPUT_TSV"
echo "GFF3 Output: $OUTPUT_GFF"
echo "End: $(date)"
