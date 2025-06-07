#!/usr/bin/env bash

#SBATCH -J trim-argania
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --time=1:50:00
#SBATCH --constraint=cal
#SBATCH --error=trim-argania.%J.err
#SBATCH --output=trim-argania.%J.out

# === Module Loading ===
module load trimmomatic/0.39

# === File Paths ===
BASE_DIR="/mnt/home/users/pab_001_uma/nhabib/Argania-spinosa-sp/Argania_spinosa"

declare -A SAMPLES
SAMPLES["1"]="ARGAN-2-D1-N1201/ARGAN-2-D1-N1201_S10"
SAMPLES["2"]="ARGAN-3-D1-N1202/ARGAN-3-D1-N1202_S11"
SAMPLES["3"]="ARGAN-5-N1203/ARGAN-5-N1203_S12"

# === Output Directory ===
OUT_DIR="trimmed_reads"
mkdir -p "$OUT_DIR"

# === Trimmomatic Parameters ===
THREADS=6
PHRED="33"
TRIM_OPTIONS="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:30 MINLEN:36"

# === Main Loop ===
for SAMPLE_NUM in "${!SAMPLES[@]}"; do
    RAW_PATH="${SAMPLES[$SAMPLE_NUM]}"
    FORWARD="$BASE_DIR/$RAW_PATH_R1_001.fastq.gz"
    REVERSE="$BASE_DIR/$RAW_PATH_R2_001.fastq.gz"

    F_OUT_P="$OUT_DIR/sample${SAMPLE_NUM}_R1_001_paired.fastq.gz"
    F_OUT_U="$OUT_DIR/sample${SAMPLE_NUM}_R1_001_unpaired.fastq.gz"
    R_OUT_P="$OUT_DIR/sample${SAMPLE_NUM}_R2_001_paired.fastq.gz"
    R_OUT_U="$OUT_DIR/sample${SAMPLE_NUM}_R2_001_unpaired.fastq.gz"

    echo "▶️ Trimming Sample ${SAMPLE_NUM}:"
    trimmomatic PE -threads $THREADS -phred$PHRED \
        "$FORWARD" "$REVERSE" \
        "$F_OUT_P" "$F_OUT_U" "$R_OUT_P" "$R_OUT_U" \
        $TRIM_OPTIONS

    echo "✅ Sample ${SAMPLE_NUM} done."
done

echo "All samples processed."
