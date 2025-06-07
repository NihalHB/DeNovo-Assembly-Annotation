#!/usr/bin/env bash
#SBATCH -J argania-qc                  # Job name
#SBATCH --ntasks=1                     # Number of tasks (processes)
#SBATCH --cpus-per-task=8             # CPUs per task
#SBATCH --mem=16gb                     # Memory per node
#SBATCH --time=1:50:00                 # Time limit
#SBATCH --constraint=cal              # Node constraint (optional)
#SBATCH --error=logs/qc.%J.err         # Standard error log
#SBATCH --output=logs/qc.%J.out        # Standard output log

# ----------------------------
# Load required modules
# ----------------------------
module load fastqc/0.11.9
module load multiqc/1.13a

# ----------------------------
# Define variables
# ----------------------------
RAW_DATA_DIR="/mnt/home/users/pab_001_uma/nhabib/Argania-spinosa-sp/Argania_spinosa"
OUT_DIR="/mnt/home/users/pab_001_uma/nhabib/Argania-spinosa-sp/01_fastqc"
THREADS=8

# ----------------------------
# Create output and log directories
# ----------------------------
mkdir -p "$OUT_DIR"
mkdir -p logs

# ----------------------------
# Run FastQC on all FASTQ files
# ----------------------------
fastqc "${RAW_DATA_DIR}"/ARGAN-2-D1-N1201/*.fastq.gz \
       "${RAW_DATA_DIR}"/ARGAN-3-D1-N1202/*.fastq.gz \
       "${RAW_DATA_DIR}"/ARGAN-5-N1203/*.fastq.gz \
       -o "$OUT_DIR" -t "$THREADS"

# ----------------------------
# Run MultiQC to summarize
# ----------------------------
multiqc "$OUT_DIR" -o "$OUT_DIR"
