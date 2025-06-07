#!/usr/bin/env bash
#SBATCH -J rnaspades-assembly
#SBATCH --ntasks=32
#SBATCH --mem=128gb
#SBATCH --time=96:00:00
#SBATCH --constraint=cal
#SBATCH --error=logs/spades.%J.err
#SBATCH --output=logs/spades.%J.out

# ----------------------------
# Load modules
# ----------------------------
module purge
module load python/3.11.4
module load spades/3.15.4

# ----------------------------
# Arguments (or use defaults)
# ----------------------------
INPUT_DIR="${1:-/mnt/home/users/pab_001_uma/nhabib/Argania-spinosa-sp/trimmomatic-argania/trimmed_reads}"
OUTPUT_DIR="${2:-/mnt/home/users/pab_001_uma/nhabib/Argania-spinosa-sp/spades-output}"
THREADS="${3:-32}"
MEMORY="${4:-128}"

mkdir -p "$OUTPUT_DIR"
mkdir -p logs

# ----------------------------
# Auto-detect input pairs
# ----------------------------
FASTQ_R1=($(find "$INPUT_DIR" -type f -name "*_R1*_paired.fastq.gz" | sort))
FASTQ_R2=($(find "$INPUT_DIR" -type f -name "*_R2*_paired.fastq.gz" | sort))

# Check if input files were found
if [ ${#FASTQ_R1[@]} -eq 0 ] || [ ${#FASTQ_R1[@]} -ne ${#FASTQ_R2[@]} ]; then
  echo "Error: Mismatched or missing paired-end files." >&2
  exit 1
fi

# ----------------------------
# Build SPAdes input options
# ----------------------------
SPADES_INPUT=""
for i in "${!FASTQ_R1[@]}"; do
  idx=$((i+1))
  SPADES_INPUT+="--pe${idx}-1 ${FASTQ_R1[$i]} --pe${idx}-2 ${FASTQ_R2[$i]} "
done

# ----------------------------
# Run rnaSPAdes
# ----------------------------
LOGFILE="$OUTPUT_DIR/rnaspades.log"

eval rnaspades.py \
  $SPADES_INPUT \
  -o "$OUTPUT_DIR" \
  --threads "$THREADS" \
  --memory "$MEMORY" \
  --disable-gzip-output \
  --ss rf \
  --phred-offset 33 \
  --debug \
  &> "$LOGFILE"
