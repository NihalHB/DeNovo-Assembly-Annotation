# Load libraries
library(Biostrings)
library(ggplot2)
library(ggpubr)
library(readr)
library(dplyr)

# ---- FILE PATHS ----
pep_file <- "transcripts_tpm_filtered.fasta.transdecoder.pep"
cds_file <- "transcripts_tpm_filtered.fasta.transdecoder.cds"
bed_file <- "transcripts_tpm_filtered.fasta.transdecoder.bed"
gff_file <- "transcripts_tpm_filtered.fasta.transdecoder.gff3"
output_pdf <- "transdecoder_report_extended.pdf"

# ---- 1. ORF AA Lengths ----
if (file.exists(pep_file)) {
  aa_seqs <- readAAStringSet(pep_file)
  orf_lengths <- width(aa_seqs)
  orf_df <- data.frame(length = orf_lengths)
  
  plot_orf_length <- ggplot(orf_df, aes(x = length)) +
    geom_histogram(bins = 50, fill = "#2c7bb6", color = "black") +
    geom_vline(aes(xintercept = median(length)), color = "red", linetype = "dashed") +
    labs(title = "ORF Lengths (AA)", x = "Length", y = "Count") +
    theme_minimal()
  
  cat("✔ ORF AA stats — min:", min(orf_lengths), " max:", max(orf_lengths), 
      " median:", median(orf_lengths), "\n")
} else {
  warning("PEP file not found.")
}

# ---- 2. CDS Nucleotide Lengths ----
if (file.exists(cds_file)) {
  cds_seqs <- readDNAStringSet(cds_file)
  cds_lengths <- width(cds_seqs)
  cds_df <- data.frame(length = cds_lengths)
  
  plot_cds_length <- ggplot(cds_df, aes(x = length)) +
    geom_histogram(bins = 50, fill = "#fdae61", color = "black") +
    geom_vline(aes(xintercept = median(length)), color = "red", linetype = "dashed") +
    labs(title = "CDS Lengths (bp)", x = "Length", y = "Count") +
    theme_minimal()
  
  cat("✔ CDS stats — min:", min(cds_lengths), " max:", max(cds_lengths), 
      " median:", median(cds_lengths), "\n")
} else {
  warning("CDS file not found.")
}

# ---- 3. BED: Transcript Position Summary ----
if (file.exists(bed_file)) {
  bed_lines <- readLines(bed_file)
  
  # Split by tab or spaces
  bed_split <- strsplit(bed_lines, "[\t ]+")
  line_lengths <- sapply(bed_split, length)
  cat("🧪 BED field counts:\n")
  print(table(line_lengths))
  
  # Use 12-column BED (BED12 format)
  valid_bed <- bed_split[line_lengths == 12]
  
  if (length(valid_bed) > 0) {
    bed_df <- as.data.frame(do.call(rbind, valid_bed), stringsAsFactors = FALSE)
    colnames(bed_df) <- c("transcript", "start", "end", "orf_id", "score", "strand",
                          "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")
    
    bed_df$start <- as.numeric(bed_df$start)
    bed_df$end <- as.numeric(bed_df$end)
    bed_df$length <- bed_df$end - bed_df$start + 1
    
    plot_bed_length <- ggplot(bed_df, aes(x = length)) +
      geom_histogram(bins = 50, fill = "#abd9e9", color = "black") +
      geom_vline(aes(xintercept = median(length)), color = "red", linetype = "dashed") +
      labs(title = "ORF Genomic Range (BED12)", x = "Length (bp)", y = "Count") +
      theme_minimal()
    
    plot_bed_strand <- ggplot(bed_df, aes(x = strand)) +
      geom_bar(fill = "#d7191c", color = "black") +
      labs(title = "Strand Usage (BED12)", x = "Strand", y = "Count") +
      theme_minimal()
    
    cat("✔ BED ORF range — min:", min(bed_df$length), 
        " max:", max(bed_df$length), 
        " median:", median(bed_df$length), "\n")
    
  } else {
    warning("No valid 12-column lines found in BED file.")
  }
  
} else {
  warning("BED file not found.")
}


# ---- EXPORT ALL TO PDF ----
pdf(output_pdf, width = 10, height = 8)
ggarrange(
  if (exists("plot_orf_length")) plot_orf_length else NULL,
  if (exists("plot_cds_length")) plot_cds_length else NULL,
  if (exists("plot_bed_length")) plot_bed_length else NULL,
  if (exists("plot_bed_strand")) plot_bed_strand else NULL,
  ncol = 2, nrow = 2
)
dev.off()

cat("Report saved to:", output_pdf, "\n")
