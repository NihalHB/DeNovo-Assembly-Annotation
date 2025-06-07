#!/usr/bin/env bash
#SBATCH -J blastx-argania
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --time=39:00:00
#SBATCH --constraint=cal
#SBATCH --output=logs/blastx.%J.out
#SBATCH --error=logs/blastx.%J.err

# ----------------------------
# Load Modules
# ----------------------------
module purge
module load blast_plus/2.15.0+

# ----------------------------
# Input Arguments or Defaults
# ----------------------------
QUERY="${1:-/mnt/home/users/pab_001_uma/nhabib/Argania-spinosa-sp/cdhit-multi/transcripts_c0.95_longest_per_cluster.fasta}"
DB="${2:-/mnt/home/users/pab_001_uma/nhabib/Argania-spinosa-sp/swissprot/swissprot}"
OUTDIR="${3:-/mnt/home/users/pab_001_uma/nhabib/Argania-spinosa-sp/blastx}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "$OUTDIR"
mkdir -p logs

OUT="$OUTDIR/blastx_output.tsv"

# ----------------------------
# Logging Info
# ----------------------------
echo "Running BLASTX at $(date)"
echo "Query:  $QUERY"
echo "DB:     $DB"
echo "Output: $OUT"
echo "Threads: $THREADS"

# ----------------------------
# Run BLASTX
# ----------------------------
blastx \
  -query "$QUERY" \
  -db "$DB" \
  -out "$OUT" \
  -evalue 1e-5 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -num_threads "$THREADS"

# ----------------------------
# Done
# ----------------------------
echo "BLASTX completed at $(date)"
echo "Output saved to: $OUT"
