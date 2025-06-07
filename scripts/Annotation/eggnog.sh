#!/usr/bin/env bash
#SBATCH -J eggnog-argania
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=72:00:00
#SBATCH --constraint=cal
#SBATCH --output=logs/eggnogmapper.%J.out
#SBATCH --error=logs/eggnogmapper.%J.err

# -------------------------------
# Load Modules
# -------------------------------
module purge
module load eggnog-mapper/2.1.12
module load diamond/2.1.8

set -e
trap 'echo "EggNOG-mapper failed at line $LINENO"; exit 1' ERR

# -------------------------------
# Input Arguments or Defaults
# -------------------------------
INPUT="${1:-/mnt/home/users/pab_001_uma/nhabib/Argania-spinosa-sp/transdecoder-output/transcripts_tpm_filtered.cleaned.pep}"
OUTDIR="${2:-/mnt2/fscratch/users/pab_001_uma/nhabib/eggnog-results}"
EGGNOG_DB="${3:-/mnt/home/users/pab_001_uma/nhabib/Argania-spinosa-sp/eggnog_data}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$OUTDIR" logs

BASENAME=$(basename "$INPUT")
PREFIX="${BASENAME%.*}"

# -------------------------------
# Run EggNOG-mapper
# -------------------------------
echo "🔍 Running EggNOG-mapper..."
echo "Input:     $INPUT"
echo "Output:    $OUTDIR"
echo "DB:        $EGGNOG_DB"
echo "Threads:   $THREADS"
echo "Started:   $(date)"

emapper.py \
  -i "$INPUT" \
  --itype proteins \
  --cpu "$THREADS" \
  -m diamond \
  --output "${PREFIX}_eggnog" \
  --output_dir "$OUTDIR" \
  --data_dir "$EGGNOG_DB"

echo "EggNOG-mapper completed at $(date)"
echo "Results saved to: $OUTDIR"
