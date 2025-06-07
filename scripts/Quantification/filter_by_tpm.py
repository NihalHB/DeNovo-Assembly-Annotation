#!/usr/bin/env python3

import os
import argparse
from collections import defaultdict
from Bio import SeqIO

# ----------------------------
# Argument Parser
# ----------------------------
parser = argparse.ArgumentParser(
    description="Filter transcripts with mean TPM ≥ threshold from Salmon quant.sf files and extract from a FASTA file."
)
parser.add_argument("-q", "--quant_dir", required=True,
                    help="Directory with subfolders containing quant.sf files")
parser.add_argument("-f", "--fasta", required=True,
                    help="Input transcriptome FASTA file")
parser.add_argument("-o", "--output", default="transcripts_tpm_filtered.fasta",
                    help="Output FASTA file name")
parser.add_argument("-t", "--threshold", type=float, default=1.0,
                    help="Mean TPM threshold (default: 1.0)")

args = parser.parse_args()

quant_dir = args.quant_dir
fasta_file = args.fasta
output_fasta = args.output
threshold = args.threshold

id_list_file = output_fasta.replace(".fasta", "_retained_ids.txt")
summary_file = output_fasta.replace(".fasta", "_summary.txt")

# ----------------------------
# Step 1: Gather TPMs
# ----------------------------
print(f"🔍 Reading TPMs from: {quant_dir}")
tpm_data = defaultdict(list)
samples_found = 0

for sample in os.listdir(quant_dir):
    sample_path = os.path.join(quant_dir, sample, "quant.sf")
    if os.path.isfile(sample_path):
        samples_found += 1
        with open(sample_path) as f:
            header = f.readline().strip().split("\t")
            name_idx = header.index("Name")
            tpm_idx = header.index("TPM")
            for line in f:
                fields = line.strip().split("\t")
                tid = fields[name_idx]
                tpm = float(fields[tpm_idx])
                tpm_data[tid].append(tpm)

print(f"Found TPMs for {len(tpm_data)} transcripts across {samples_found} samples")

# ----------------------------
# Step 2: Filter transcripts
# ----------------------------
filtered_ids = [tid for tid, tpms in tpm_data.items() if (sum(tpms) / len(tpms)) >= threshold]
print(f"Retaining {len(filtered_ids)} transcripts with mean TPM ≥ {threshold}")

# ----------------------------
# Step 3: Extract sequences
# ----------------------------
print(f"🧬 Extracting from: {fasta_file}")
seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
written = 0

with open(output_fasta, "w") as out_fa, open(id_list_file, "w") as idf:
    for tid in filtered_ids:
        if tid in seq_dict:
            SeqIO.write(seq_dict[tid], out_fa, "fasta")
            idf.write(tid + "\n")
            written += 1
        else:
            print(f"Warning: {tid} not found in FASTA")

print(f"{written} sequences written to {output_fasta}")
print(f"Retained transcript IDs saved to {id_list_file}")

# ----------------------------
# Step 4: Write Summary
# ----------------------------
total_input = len(seq_dict)
percent_retained = (written / total_input) * 100

with open(summary_file, "w") as sf:
    sf.write("=== TPM Filtering Summary ===\n")
    sf.write(f"Total input transcripts:        {total_input}\n")
    sf.write(f"TPM threshold:                  {threshold}\n")
    sf.write(f"Transcripts retained (TPM ≥ t): {len(filtered_ids)}\n")
    sf.write(f"Transcripts written to FASTA:   {written}\n")
    sf.write(f"Transcripts dropped:            {total_input - written}\n")
    sf.write(f"Percent retained:               {percent_retained:.2f}%\n")

print(f"📊 Summary written to: {summary_file}")
