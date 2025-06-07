**Comprehensive Pipeline for De Novo Assembly and Functional Annotation of non-model organisms: case study on the Argania spinosa Transcriptome**

---

## Project Overview

This repository provides a complete and reproducible pipeline for de novo transcriptome analysis of *Argania spinosa*, a non-model organism that is an ecologically and economically important tree species endemic to southwestern Morocco. The pipeline includes quality control, assembly, redundancy reduction, open reading frame prediction, functional annotation, repeat detection, transcript quantification, and machine learning-based classification.

This workflow integrates Bash, R, and Python scripts for each module and is optimized for execution on HPC environments or local UNIX-based systems.

---

## 🔄 Pipeline Workflow Structure

```
DeNovo-Assembly-Annotation/
├── scripts/
│   ├── QC/                     # FastQC, Trimmomatic
│   ├── Assembly/               # RNA-SPAdes, CD-HIT
│   ├── orf-prediction/         # TransDecoder
│   ├── Annotation/             # BLASTx, InterProScan, RepeatMasker
│   ├── Quantification/         # TPM filtering
│   ├── ML/                     # Machine learning
│   └── Visualizations/         # Plotting scripts
|
└── run_all.sh                 # Master automation script
```

---

## How to Clone and Set Up

```bash
git clone https://github.com/NihalHB/DeNovo-Assembly-Annotation.git
cd DeNovo-Assembly-Annotation
conda env create -f environment.yml
conda activate argania-pipeline
```


## How to Run Entire Pipeline

```bash
bash run_all.sh
```

This will execute the following steps:

1. Quality Control (FastQC, Trimmomatic)
2. Assembly (rnaSPAdes)
3. Clustering (CD-HIT)
4. ORF Prediction (TransDecoder)
5. Annotation (BLASTx, InterProScan, RepeatMasker, EggNOG)
6. TPM Filtering
7. Machine Learning Classification
8. Visualizations

---

Make executable:

```bash
chmod +x run_all.sh
```

---

## 📊 Outputs

* **Trimmed Reads**
* **Assembled Transcripts**
* **Clustered Contigs**
* **Predicted ORFs**
* **Annotations**
* **TPM Filtered Files**
* **ML Classification Results**
* **Plots and Visual Reports**

---

## Contact

Nihal Habib
[nhabib@um6ss.ma](mailto:nhabib@um6ss.ma)
[@NihalHB](https://github.com/NihalHB)

---

## License

MIT License
