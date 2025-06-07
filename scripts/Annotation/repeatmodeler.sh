#!/usr/bin/env bash
#SBATCH -J repeatmodeler-argania
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=72:00:00
#SBATCH --output=logs/repeatmodeler.%J.out
#SBATCH --error=logs/repeatmodeler.%J.err
#SBATCH --constraint=cal

# -------------------------------
# Load Modules
# -------------------------------
module purge
module load repeatmodeler/2.0.4_new
module load repeatmasker/4.1.5

set -e
trap 'echo "Error at line $LINENO. Exiting."; exit 1' ERR

# -------------------------------
# Define Paths
# -------------------------------
BASE="/mnt/home/users/pab_001_uma/nhabib/Argania-spinosa-sp"
INPUT="${1:-$BASE/transdecoder-output/transcripts_tpm_filtered.fasta}"
OUTDIR="${2:-$BASE/repeatmodeler_output}"
DBNAME="${3:-argania_db}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "$OUTDIR" logs
cd "$OUTDIR"

# -------------------------------
# Build Database
# -------------------------------
echo "📦 Building RepeatModeler database at $(date)"
BuildDatabase -name "$DBNAME" -engine ncbi "$INPUT"

# -------------------------------
# Run RepeatModeler
# -------------------------------
echo "🧬 Running RepeatModeler with $THREADS threads..."
RepeatModeler -database "$DBNAME" -engine ncbi -threads "$THREADS"

echo "RepeatModeler completed at $(date)"
echo "Results located in: $OUTDIR"
