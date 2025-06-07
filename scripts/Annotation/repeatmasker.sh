#!/usr/bin/env bash
#SBATCH -J repeatmasker-argania
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=72:00:00
#SBATCH --constraint=cal
#SBATCH --output=logs/repeatmasker.%J.out
#SBATCH --error=logs/repeatmasker.%J.err

# -------------------------------
# Load Module
# -------------------------------
module purge
module load repeatmasker/4.1.5

set -e
trap 'echo "Error at line $LINENO. Exiting."; exit 1' ERR

# -------------------------------
# Input Arguments or Defaults
# -------------------------------
BASE="${1:-/mnt2/fscratch/users/pab_001_uma/nhabib}"
INPUT="${2:-$BASE/transcripts_cdhit.fasta}"
LIBRARY="${3:-$BASE/RM_27332.WedJun4162505202/consensi.fa.classified}"
OUTDIR="${4:-$BASE/Argania-annotation-outputs/repeatmasker-results}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "$OUTDIR" logs

# -------------------------------
# Run RepeatMasker
# -------------------------------
echo "Running RepeatMasker..."
echo "Input:    $INPUT"
echo "Library:  $LIBRARY"
echo "Output:   $OUTDIR"
echo "Threads:  $THREADS"
echo "Started:  $(date)"

RepeatMasker \
  -pa "$THREADS" \
  -lib "$LIBRARY" \
  -gff \
  -dir "$OUTDIR" \
  "$INPUT"

echo "RepeatMasker completed at: $(date)"
echo "Output directory: $OUTDIR"
