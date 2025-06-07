#!/usr/bin/env bash
#SBATCH -J cdhit-argania-multi
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=48:00:00
#SBATCH --constraint=cal
#SBATCH --output=logs/cdhit-multi.%J.out
#SBATCH --error=logs/cdhit-multi.%J.err

# ----------------------------
# Load CD-HIT
# ----------------------------
module purge
module load cdhit/4.8.1

# ----------------------------
# Input/output arguments or defaults
# ----------------------------
INPUT="${1:-/mnt/home/users/pab_001_uma/nhabib/Argania-spinosa-sp/spades-output/transcripts.fasta}"
OUTDIR="${2:-/mnt/home/users/pab_001_uma/nhabib/Argania-spinosa-sp/cdhit-multi}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
THRESHOLDS=(0.98 0.95 0.90 0.85)

mkdir -p "$OUTDIR"
mkdir -p logs

# ----------------------------
# Header for summary
# ----------------------------
SUMMARY_FILE="$OUTDIR/cluster_summary.csv"
echo "Threshold,Num_Clusters" > "$SUMMARY_FILE"

# ----------------------------
# Run CD-HIT-EST for each threshold
# ----------------------------
for C in "${THRESHOLDS[@]}"; do
    echo ">>> Running CD-HIT-EST at -c $C"
    
    OUT="${OUTDIR}/transcripts_c${C}"
    
    cd-hit-est \
      -i "$INPUT" \
      -o "${OUT}.fasta" \
      -c "$C" \
      -n 10 \
      -G 1 \
      -aS 0.8 \
      -aL 0.8 \
      -d 0 \
      -M 64000 \
      -T "$THREADS"
    
    # Count clusters
    CLSTR="${OUT}.fasta.clstr"
    if [[ -f "$CLSTR" ]]; then
      NUM_CLUSTERS=$(grep -c "^>Cluster" "$CLSTR")
      echo ">>> Threshold $C → $NUM_CLUSTERS clusters"
      echo "$C,$NUM_CLUSTERS" >> "$SUMMARY_FILE"
    else
      echo "⚠️  Warning: $CLSTR not found. Skipping summary for -c $C"
    fi
done

echo "CD-HIT runs complete. Summary saved to: $SUMMARY_FILE"
