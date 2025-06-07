# -------------------------------
# Load Required Packages
# -------------------------------
packages <- c("tidyverse", "ggplot2", "scales")
lapply(packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)
})

# -------------------------------
# Define File Paths
# -------------------------------
blast_file <- "blastx_output.tsv"
lengths_file <- "transcript_lengths.csv"
output_dir <- "/Users/Apple/Desktop/projects/New argania analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------------------
# Load BLAST and Length Data
# -------------------------------
col_names <- c("qseqid", "sseqid", "pident", "aln_length", "mismatch", "gapopen",
               "qstart", "qend", "sstart", "send", "evalue", "bitscore")

blast <- read_tsv(blast_file, col_names = col_names, show_col_types = FALSE)
lengths <- read_csv(lengths_file, show_col_types = FALSE)

# -------------------------------
# Calculate Query Coverage (%)
# -------------------------------
blast_cov <- blast %>%
  left_join(lengths, by = c("qseqid" = "Transcript")) %>%
  mutate(
    coverage = round(((qend - qstart + 1) / Length) * 100, 2)
  ) %>%
  filter(!is.na(coverage))  # remove rows with missing length info

# -------------------------------
# Filter Valid Coverage (0–100%)
# -------------------------------
blast_clean <- blast_cov %>%
  filter(coverage >= 0 & coverage <= 100)

# Save cleaned table
write.csv(blast_clean, file.path(output_dir, "blastx_cleaned_with_coverage.csv"), row.names = FALSE)

# -------------------------------
# Identify Top Fully Covered Transcripts
# -------------------------------
top_full_hits <- blast_clean %>%
  filter(coverage == 100) %>%
  arrange(desc(bitscore))

write.csv(top_full_hits, file.path(output_dir, "top_full_coverage_hits.csv"), row.names = FALSE)

# -------------------------------
# Plot 1: Histogram of Cleaned Coverage
# -------------------------------
p1 <- ggplot(blast_clean, aes(x = coverage)) +
  geom_histogram(fill = "#66c2a5", bins = 50, color = "black") +
  geom_vline(aes(xintercept = median(coverage, na.rm = TRUE)), color = "red", linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(title = "Histogram of Transcript Query Coverage", x = "Coverage (%)", y = "Frequency")

ggsave(file.path(output_dir, "blastx_coverage_histogram.png"), p1, dpi = 300, width = 8, height = 5)

# -------------------------------
# Plot 2: Violin Plot by pident Bin
# -------------------------------
blast_clean$pident_bin <- cut(blast_clean$pident,
                              breaks = c(0, 50, 70, 85, 95, 100),
                              labels = c("<50%", "50–70%", "70–85%", "85–95%", "95–100%"),
                              include.lowest = TRUE)

p2 <- ggplot(blast_clean, aes(x = pident_bin, y = coverage)) +
  geom_violin(fill = "#2c7bb6", color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Coverage by % Identity", x = "% Identity Bin", y = "Coverage (%)")

ggsave(file.path(output_dir, "blastx_coverage_violin_plot.png"), p2, dpi = 300, width = 7, height = 5)

# -------------------------------
# Plot 3: Density Plot of Coverage
# -------------------------------
p3 <- ggplot(blast_clean, aes(x = coverage)) +
  geom_density(fill = "#fdae61", alpha = 0.6) +
  theme_minimal(base_size = 14) +
  labs(title = "Density of Transcript Query Coverage", x = "Coverage (%)", y = "Density")

ggsave(file.path(output_dir, "blastx_coverage_density_plot.png"), p3, dpi = 300, width = 7, height = 5)


# -------------------------------
# Extract UniProt ID, Protein Name, and Species
# -------------------------------
blast <- blast %>%
  mutate(
    uniprot_id = str_extract(subject_id, "sp\\|([^|]+)\\|") %>% str_remove_all("sp\\||\\|"),
    protein_name = str_extract(subject_id, "\\|[^|]+\\|([^_]+)_") %>% str_remove_all("\\|"),
    species = str_extract(subject_id, "_([A-Z0-9]+)$") %>% str_remove("_")
  )

# View parsed annotation table (remove duplicates if needed)
annotation_table <- blast %>%
  select(query_id, uniprot_id, protein_name, species, identity, evalue, bit_score) %>%
  distinct()

# ---- 4. Count Hits per Gene/Protein Family ----
gene_family_counts <- annotation_table %>%
  count(protein_name, sort = TRUE)

# -------------------------------
# Plot 4: Top 10 most common protein families
# -------------------------------
ggplot(head(gene_family_counts, 50), aes(x = reorder(protein_name, n), y = n)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 50 Protein Families in BLASTx Hits", x = "Protein Family", y = "Number of Hits") +
  theme_minimal()

# ----  Species Distribution ----
species_counts <- annotation_table %>%
  count(species, sort = TRUE)

# -------------------------------
# Plot 5: Top 10 species
# -------------------------------
ggplot(head(species_counts, 50), aes(x = reorder(species, n), y = n)) +
  geom_col(fill = "darkgreen") +
  coord_flip() +
  labs(title = "Top 50 Species in BLASTx Hits", x = "Species", y = "Number of Hits") +
  theme_minimal()

# ---- Save Results ----
write.csv(annotation_table, "blastx_annotation_table.csv", row.names = FALSE)
write.csv(gene_family_counts, "blastx_gene_family_counts.csv", row.names = FALSE)
write.csv(species_counts, "blastx_species_counts.csv", row.names = FALSE)

# -------------------------------
# Done
# -------------------------------
cat("All outputs saved to:", output_dir, "\n")
