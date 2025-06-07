#!/usr/bin/env bash
#SBATCH -J salmon-quant
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --time=39:00:00
#SBATCH --constraint=cal
#SBATCH --output=logs/salmon_quant.%J.out
#SBATCH --error=logs/salmon_quant.%J.err

# ----------------------------
# Load module
# ----------------------------
module purge
module load salmon/1.10.2

# ----------------------------
# Define paths and arguments
# ----------------------------
BASE_DIR="${1:-/mnt/home/users/pab_001_uma/nhabib/Argania-spinosa-sp}"
READS_DIR="${BASE_DIR}/Argania_spinosa"
SALMON_INDEX="${BASE_DIR}/salmon_index"
OUTPUT_DIR="${BASE_DIR}/quantification"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# ----------------------------
# Loop through each sample
# ----------------------------
for SAMPLE_DIR in "${READS_DIR}"/ARGAN-*; do
    SAMPLE_NAME=$(basename "$SAMPLE_DIR")

    R1=$(find "$SAMPLE_DIR" -type f -name "*_R1_001.fastq.gz")
    R2=$(find "$SAMPLE_DIR" -type f -name "*_R2_001.fastq.gz")

    if [[ -z "$R1" || -z "$R2" ]]; then
        echo "Skipping $SAMPLE_NAME — missing R1 or R2 files"
        continue
    fi

    SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE_NAME}"
    mkdir -p "$SAMPLE_OUT"

    echo "🔄 Running Salmon on $SAMPLE_NAME..."

    salmon quant -i "$SALMON_INDEX" -l A \
        -1 "$R1" -2 "$R2" \
        -p "$THREADS" \
        -o "$SAMPLE_OUT"

    echo "Completed: $SAMPLE_NAME → $SAMPLE_OUT"
done

echo "All samples processed!"
