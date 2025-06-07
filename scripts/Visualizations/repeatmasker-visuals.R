# ==========================
# RepeatMasker Output Parser and Visualizer
# ==========================

library(tidyverse)

# ---- Define Paths ----
repeat_dir <- "."  # Change if needed
output_dir <- file.path(repeat_dir, "repeatmasker_plots")
dir.create(output_dir, showWarnings = FALSE)

rm_out_file <- file.path(repeat_dir, "transcripts_tpm_filtered.fasta.out")

# ---- Read and Clean ----
rm_lines <- readLines(rm_out_file)
rm_hits <- rm_lines[grep("^\\s*\\d+", rm_lines)]

# Split using regex: multiple spaces or tabs
rm_df <- rm_hits %>%
  str_squish() %>%
  str_split(" ", simplify = FALSE) %>%
  map_dfr(~ tibble(
    SW_score = .x[1],
    perc_div = .x[2],
    perc_del = .x[3],
    perc_ins = .x[4],
    query_seq = .x[5],
    query_start = .x[6],
    query_end = .x[7],
    query_left = .x[8],
    strand = .x[9],
    repeat_name = .x[10],
    repeat_class_family = .x[11],
    repeat_start = .x[12],
    repeat_end = .x[13],
    repeat_left = .x[14],
    ID = .x[15]
  ))

# ---- Clean and Parse Columns ----
rm_df <- rm_df %>%
  mutate(across(everything(), ~ str_replace_all(.x, "[()]", ""))) %>%
  mutate(across(c(SW_score, perc_div, perc_del, perc_ins,
                  query_start, query_end, repeat_start, repeat_end, ID),
                as.numeric)) %>%
  mutate(repeat_length = abs(query_end - query_start) + 1) %>%
  separate(repeat_class_family, into = c("Class", "Family"), sep = "/", fill = "right")

# ---- Save Cleaned CSV ----
write_csv(rm_df, file.path(output_dir, "Parsed_RepeatMasker_Output.csv"))

# ---- Pie Chart: Repeat Class Composition ----
class_summary <- rm_df %>%
  count(Class, sort = TRUE) %>%
  filter(!is.na(Class))

ggplot(class_summary, aes(x = "", y = n, fill = Class)) +
  geom_col(width = 1, color = "white") +
  coord_polar("y") +
  labs(title = "Repeat Class Composition", fill = "Class") +
  theme_void() +
  ggsave(file.path(output_dir, "Repeat_Classes_PieChart.png"), dpi = 300, width = 6, height = 6)

# ---- Bar Plot: Top 10 Repeat Families ----
top_families <- rm_df %>%
  count(Family, sort = TRUE) %>%
  slice_max(n, n = 20) %>%
  filter(!is.na(Family))

ggplot(top_families, aes(x = reorder(Family, n), y = n)) +
  geom_col(fill = "#2c7bb6") +
  coord_flip() +
  labs(title = "Top 20 Repeat Families", x = "Family", y = "Count") +
  theme_minimal() +
  ggsave(file.path(output_dir, "Top20_Repeat_Families_Barplot.png"), dpi = 300, width = 7, height = 5)

# ---- Bar Plot: Repeat Length Contribution by Class ----
length_summary <- rm_df %>%
  filter(!is.na(Class)) %>%
  group_by(Class) %>%
  summarise(total_length = sum(repeat_length, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_length))

ggplot(length_summary, aes(x = reorder(Class, total_length), y = total_length)) +
  geom_col(fill = "#d7191c") +
  coord_flip() +
  labs(title = "Total Length of Repeats by Class", x = "Repeat Class", y = "Total Repeat Length (bp)") +
  theme_minimal() +
  ggsave(file.path(output_dir, "Repeat_Length_By_Class.png"), dpi = 300, width = 7, height = 5)

cat("✅ RepeatMasker output parsed successfully into:", output_dir, "\n")

# Load libraries
library(tidyverse)

# Read parsed RepeatMasker output
rm_df <- read_csv("repeatmasker_plots/Parsed_RepeatMasker_Output.csv")

# Filter to Simple Repeats
simple_df <- rm_df %>%
  filter(Class == "Simple_repeat") %>%
  mutate(
    motif = str_remove_all(repeat_name, "[()]"),  # remove parentheses
    
    repeat_type = case_when(
      str_length(motif) == 1 ~ "Mono",
      str_length(motif) == 2 ~ "Di",
      str_length(motif) == 3 ~ "Tri",
      str_length(motif) == 4 ~ "Tetra",
      str_length(motif) > 4 ~ "Rich",
      TRUE ~ "Other"
    )
  )

# Optional: keep only most common motifs to avoid clutter
top_motifs <- simple_df %>%
  count(motif, sort = TRUE) %>%
  slice_head(n = 15) %>%  # adjust number to show top 10–20 motifs
  pull(motif)

filtered_df <- simple_df %>%
  filter(motif %in% top_motifs)

# Count for stacked bar
repeat_counts <- filtered_df %>%
  count(repeat_type, motif)

# Plot
ggplot(repeat_counts, aes(x = repeat_type, y = n, fill = motif)) +
  geom_col(position = "stack") +
  labs(
    title = "Stacked Bar Plot of Simple Repeat Types (from RepeatMasker)",
    x = "Repeat Type", y = "Count", fill = "Motif"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

# Save plot
ggsave("RepeatMasker_Simple_Repeat_Stacked_TopMotifs.png", dpi = 300, width = 8, height = 5)

ggplot(rm_df %>% filter(!is.na(Class)), aes(x = Class, y = repeat_length)) +
  geom_boxplot(outlier.shape = NA) +
  coord_flip() +
  labs(title = "Repeat Length Distribution by Class",
       y = "Length (bp)", x = "Repeat Class") +
  theme_minimal()
ggsave("Repeat Length Distribution by Class.png", dpi = 300, width = 8, height = 5)

repeat_density <- rm_df %>%
  count(query_seq)

ggplot(repeat_density, aes(x = n)) +
  geom_histogram(binwidth = 1, fill = "#2c7bb6") +
  labs(title = "Distribution of Repeat Counts per Transcript",
       x = "Number of Repeats", y = "Transcript Count") +
  theme_minimal()
ggsave("Distribution of Repeat Counts per Transcript.png", dpi = 300, width = 8, height = 5)

masked_summary <- rm_df %>%
  distinct(query_seq) %>%
  mutate(masked = "Repeat-Masked")

all_transcripts <- data.frame(
  query_seq = unique(rm_df$query_seq),
  masked = "Total"
)

combined <- bind_rows(masked_summary, all_transcripts) %>%
  count(masked)

ggplot(combined, aes(x = "", y = n, fill = masked)) +
  geom_col(width = 1) +
  coord_polar("y") +
  theme_void() +
  labs(title = "Proportion of Repeat-Masked Transcripts")
