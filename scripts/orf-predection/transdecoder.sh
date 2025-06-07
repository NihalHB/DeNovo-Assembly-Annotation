#!/usr/bin/env bash
#SBATCH -J transdecoder-argania
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --time=39:00:00
#SBATCH --constraint=cal
#SBATCH --output=logs/transdecoder.%J.out
#SBATCH --error=logs/transdecoder.%J.err

# -------------------------------
# Load Conda & Activate Env
# -------------------------------
module purge
module load miniconda
source /mnt/home/soft/miniconda/programs/x86_64/miniconda3_py310_23.1.0/etc/profile.d/conda.sh
conda activate transdecoder

set -e
trap 'echo "❌ Error at line $LINENO. Exiting."; exit 1' ERR

# -------------------------------
# Input/Output Setup
# -------------------------------
TRANSCRIPTS="${1:-/mnt/home/users/pab_001_uma/nhabib/Argania-spinosa-sp/quantification/transcripts_tpm_filtered.fasta}"
OUTDIR="${2:-/mnt/home/users/pab_001_uma/nhabib/Argania-spinosa-sp/transdecoder-output}"
LOGFILE="$OUTDIR/transdecoder.log"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "$OUTDIR" logs
cp "$TRANSCRIPTS" "$OUTDIR/"
cd "$OUTDIR" || exit 1

BASENAME=$(basename "$TRANSCRIPTS")

# -------------------------------
# Log Start
# -------------------------------
echo "Starting TransDecoder on $BASENAME" | tee "$LOGFILE"
start_time=$(date)

# -------------------------------
# Step 1: Long ORFs
# -------------------------------
echo "Step 1: Running TransDecoder.LongOrfs..." | tee -a "$LOGFILE"
TransDecoder.LongOrfs -t "$BASENAME" >> "$LOGFILE" 2>&1

# -------------------------------
# Step 2: Coding Prediction
# -------------------------------
echo "Step 2: Running TransDecoder.Predict..." | tee -a "$LOGFILE"
TransDecoder.Predict -t "$BASENAME" >> "$LOGFILE" 2>&1

# -------------------------------
# Done
# -------------------------------
end_time=$(date)
echo "Done. Outputs in: $OUTDIR" | tee -a "$LOGFILE"
echo "Start: $start_time" | tee -a "$LOGFILE"
echo "End:   $end_time"   | tee -a "$LOGFILE"
