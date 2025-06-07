# ============================================
# EggNOG-mapper Annotation Analysis - Argania
# ============================================

# Load required libraries
library(tidyverse)
library(ggplot2)
library(forcats)
library(scales)

# Read the EggNOG annotations file
eggnog <- read_tsv("eggnog_annotations.annotations", comment = "#", show_col_types = FALSE)

# -----------------------
# Clean and Filter Data
# -----------------------

# Remove rows without meaningful annotations
annotated <- eggnog %>%
  filter(!is.na(Description) & Description != "-")

# -----------------------
# Top Functional Descriptions
# -----------------------

top_desc <- annotated %>%
  count(Description, sort = TRUE) %>%
  slice_max(n, n = 15)

ggplot(top_desc, aes(x = fct_reorder(Description, n), y = n)) +
  geom_col(fill = "#1f78b4") +
  coord_flip() +
  labs(title = "Top Functional Descriptions",
       x = "Function", y = "Frequency") +
  theme_minimal(base_size = 12)

# -----------------------
# GO Term Frequency
# -----------------------

go_terms <- annotated %>%
  filter(!is.na(GOs)) %>%
  separate_rows(GOs, sep = ",") %>%
  count(GOs, sort = TRUE) %>%
  slice_max(n, n = 15)

ggplot(go_terms, aes(x = fct_reorder(GOs, n), y = n)) +
  geom_col(fill = "#33a02c") +
  coord_flip() +
  labs(title = "Top GO Terms", x = "GO Term", y = "Count") +
  theme_minimal(base_size = 12)

# -----------------------
# KEGG Pathway Frequency
# -----------------------

kegg <- annotated %>%
  filter(!is.na(KEGG_Pathway)) %>%
  separate_rows(KEGG_Pathway, sep = ",") %>%
  count(KEGG_Pathway, sort = TRUE) %>%
  slice_max(n, n = 15)

ggplot(kegg, aes(x = fct_reorder(KEGG_Pathway, n), y = n)) +
  geom_col(fill = "#e31a1c") +
  coord_flip() +
  labs(title = "Top KEGG Pathways", x = "KEGG Pathway", y = "Count") +
  theme_minimal(base_size = 12)

# -----------------------
# PFAM Domain Frequency
# -----------------------

pfams <- annotated %>%
  filter(!is.na(PFAMs)) %>%
  separate_rows(PFAMs, sep = ",") %>%
  count(PFAMs, sort = TRUE) %>%
  slice_max(n, n = 15)

ggplot(pfams, aes(x = fct_reorder(PFAMs, n), y = n)) +
  geom_col(fill = "#ff7f00") +
  coord_flip() +
  labs(title = "Top PFAM Domains", x = "PFAM ID", y = "Count") +
  theme_minimal(base_size = 12)

cat("✅ EggNOG-mapper annotation figures generated successfully.\n")
