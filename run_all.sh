#!/usr/bin/env bash

set -euo pipefail

echo "Step 1: Quality Control"
bash scripts/QC/trimmomatic.sh

echo "Step 2: Assembly"
bash scripts/Assembly/spades.sh

echo "Step 3: Redundancy Removal"
bash scripts/Assembly/cd-hit.sh

echo "Step 4: ORF Prediction"
bash scripts/orf-prediction/transdecoder.sh

echo "Step 5: Functional Annotation"
bash scripts/Annotation/blastx.sh
bash scripts/Annotation/interproscan.sh
bash scripts/Annotation/repeatmasker.sh
bash scripts/Annotation/eggnog.sh

echo "Step 6: Quantification"
python scripts/Quantification/filter_by_tpm.py

echo "Step 7: Machine Learning"
Rscript scripts/ML/machine-learning.R

echo "Step 8: Visualizations"
Rscript scripts/Visualizations/cdhit-visuals.R
