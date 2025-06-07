# Load required packages
library(ggplot2)
library(readr)

# Read the CSV file
summary_file <- "/mnt/home/users/pab_001_uma/nhabib/Argania-spinosa-sp/cdhit-multi/cluster_summary.csv"
cluster_data <- read_csv(summary_file, col_names = c("Identity", "NumClusters"))

# Make sure Identity is treated as numeric
cluster_data$Identity <- as.numeric(cluster_data$Identity)

# Plot
p <- ggplot(cluster_data, aes(x = Identity, y = NumClusters)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = cluster_data$Identity) +
  labs(
    title = "CD-HIT Cluster Count vs. Identity Threshold",
    x = "CD-HIT Identity Threshold (-c)",
    y = "Number of Clusters"
  ) +
  theme_minimal(base_size = 14)

# Save plot
ggsave("cdhit_cluster_summary_plot.png", p, width = 8, height = 5, dpi = 600)

# Also display
print(p)
