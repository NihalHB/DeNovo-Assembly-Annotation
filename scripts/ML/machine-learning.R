# ================================
# 1. Setup & Required Libraries
# ================================
packages <- c("tidyverse", "randomForest", "caret", "pROC", "Rtsne", "factoextra", "Biostrings", "readr", "e1071", "VennDiagram", "grid", "cowplot", "ggpubr")
installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)
if (length(to_install)) install.packages(to_install)
lapply(packages, library, character.only = TRUE)

# ================================
# 2. Input Data Preparation
# ================================
setwd("/Users/Apple/Desktop/projects/New argania analysis/")

# ORF Lengths
pep <- readAAStringSet("transcripts_tpm_filtered.fasta.transdecoder.pep")
orf_df <- data.frame(transcript = sapply(strsplit(names(pep), " "), `[`, 1), orf_length = width(pep))
write.csv(orf_df, "orf_lengths.csv", row.names = FALSE)

# Quantification
samples <- c("ARGAN-2", "ARGAN-3", "ARGAN-5")
files <- file.path("quantification", c("ARGAN-2-D1-N1201", "ARGAN-3-D1-N1202", "ARGAN-5-N1203"), "quant.sf")
names(files) <- samples
txi <- tximport(files, type = "salmon", txOut = TRUE)
tpm_df <- as.data.frame(txi$abundance)
tpm_df$transcript <- rownames(tpm_df)
tpm_df <- tpm_df[, c("transcript", samples)]
write.csv(tpm_df, "tpm_all.csv", row.names = FALSE)

# Annotation Features
repeat_gff <- read_tsv("transcripts_tpm_filtered.fasta.out.gff", comment = "#", col_names = FALSE)
te <- repeat_gff %>% filter(X3 == "dispersed_repeat") %>% select(Transcript_ID = X1) %>% distinct() %>% mutate(Has_TE = 1)
interpro <- read_tsv("ips_output_20250602_204244.tsv", col_names = FALSE, comment = "#")
domains <- interpro %>% filter(!is.na(X12)) %>% select(Transcript_ID = X1, InterPro_ID = X12) %>% distinct() %>% count(Transcript_ID, name = "Domain_Count")
orf <- read_csv("orf_lengths.csv") %>% rename(Transcript_ID = transcript) %>% mutate(Transcript_ID = str_remove(Transcript_ID, "\\.p[0-9]+$"))
tpm <- read_csv("tpm_all.csv") %>% rename(Transcript_ID = transcript) %>% pivot_longer(cols = -Transcript_ID, names_to = "Sample", values_to = "TPM") %>% group_by(Transcript_ID) %>% summarise(TPM = mean(TPM, na.rm = TRUE))
fasta <- readDNAStringSet("transcripts_tpm_filtered.fasta")
transcript_features <- tibble(
  Transcript_ID = names(fasta),
  Length = width(fasta),
  GC = rowSums(letterFrequency(fasta, c("G", "C"), as.prob = TRUE))
)

# Merge All Features
features <- transcript_features %>%
  left_join(tpm, by = "Transcript_ID") %>%
  left_join(orf, by = "Transcript_ID") %>%
  left_join(te, by = "Transcript_ID") %>%
  left_join(domains, by = "Transcript_ID") %>%
  mutate(
    Has_TE = replace_na(Has_TE, 0),
    TPM = replace_na(TPM, 0),
    ORF_Length = replace_na(orf_length, 0),
    Domain_Count = replace_na(Domain_Count, 0),
    Class = factor(if_else(Domain_Count > 0 | ORF_Length > 100, "Functional", "NonFunctional"))
  ) %>%
  drop_na(GC, Length, TPM)
write_csv(features, "argania_full_features.csv")

# ===================================
# 3. Train-Test Split (No Data Leakage)
# ===================================
set.seed(42)
train_idx <- createDataPartition(features$Class, p = 0.7, list = FALSE)
training_data <- features[train_idx, ]
test_data <- features[-train_idx, ]

# ================================
# 4. Train Classifiers
# ================================
control <- trainControl(method = "cv", number = 5, classProbs = TRUE, savePredictions = TRUE)
rf_model <- train(Class ~ Length + GC + TPM + Has_TE, data = training_data, method = "rf", trControl = control, ntree = 500)
svm_model <- train(Class ~ Length + GC + TPM + Has_TE, data = training_data, method = "svmRadial", trControl = control, preProcess = c("center", "scale"), tuneLength = 5)

# ================================
# 5. Predict & Evaluate on Holdout Test Set
# ================================
test_data$RF_Predicted <- predict(rf_model, newdata = test_data)
test_data$SVM_Predicted <- predict(svm_model, newdata = test_data)
test_data$RF_Prob_Functional <- predict(rf_model, newdata = test_data, type = "prob")$Functional
test_data$SVM_Prob_Functional <- predict(svm_model, newdata = test_data, type = "prob")$Functional

write_csv(test_data, "test_predictions.csv")

# ================================
# 6. Visualizations & PDF Summary
# ================================
# ROC Curves
rf_roc <- roc(response = test_data$Class, predictor = test_data$RF_Prob_Functional)
svm_roc <- roc(response = test_data$Class, predictor = test_data$SVM_Prob_Functional)

# Confusion Matrices
cm_rf <- confusionMatrix(test_data$RF_Predicted, test_data$Class)
cm_svm <- confusionMatrix(test_data$SVM_Predicted, test_data$Class)

# Feature Importance
importance_df <- varImp(rf_model)$importance %>% rownames_to_column("Feature")

# Metrics Table
comparison_metrics <- tibble(
  Model = c("Random Forest", "SVM"),
  Accuracy = c(cm_rf$overall["Accuracy"], cm_svm$overall["Accuracy"]),
  Kappa = c(cm_rf$overall["Kappa"], cm_svm$overall["Kappa"]),
  Sensitivity = c(cm_rf$byClass["Sensitivity"], cm_svm$byClass["Sensitivity"]),
  Specificity = c(cm_rf$byClass["Specificity"], cm_svm$byClass["Specificity"])
)
write_csv(comparison_metrics, "model_metrics_summary.csv")

# PDF Export
pdf("final_model_evaluation.pdf", width = 8, height = 10)
plot(rf_roc, col = "#d7191c", main = sprintf("ROC Curve (AUC RF: %.3f, AUC SVM: %.3f)", auc(rf_roc), auc(svm_roc)))
lines(svm_roc, col = "#2c7bb6")
legend("bottomright", legend = c("Random Forest", "SVM"), col = c("#d7191c", "#2c7bb6"), lwd = 2)

rf_table <- as.data.frame(cm_rf$table)
p1 <- ggplot(rf_table, aes(Prediction, Reference, fill = Freq)) +
  geom_tile() + geom_text(aes(label = Freq)) +
  scale_fill_gradient(low = "white", high = "#d7191c") +
  theme_minimal() + ggtitle("Confusion Matrix - Random Forest") +
  xlab("Predicted") + ylab("Actual")
print(p1)

svm_table <- as.data.frame(cm_svm$table)
p2 <- ggplot(svm_table, aes(Prediction, Reference, fill = Freq)) +
  geom_tile() + geom_text(aes(label = Freq)) +
  scale_fill_gradient(low = "white", high = "#2c7bb6") +
  theme_minimal() + ggtitle("Confusion Matrix - SVM") +
  xlab("Predicted") + ylab("Actual")
print(p2)

ggplot(importance_df, aes(x = reorder(Feature, Overall), y = Overall)) +
  geom_col(fill = "#2c7bb6") + coord_flip() + theme_minimal() +
  labs(title = "Feature Importance (Random Forest)", x = "Feature", y = "Importance")

theme_tbl <- tableGrob(round(select(comparison_metrics, -Model), 4), rows = comparison_metrics$Model)
grid.newpage(); grid.draw(theme_tbl)

dev.off()

# ================================
# 7. Precision and Recall (Functional Class)
# ================================
# Helper function to extract precision/recall for 'Functional'
get_pr <- function(cm) {
  precision <- cm$byClass["Pos Pred Value"]
  recall <- cm$byClass["Sensitivity"]
  return(c(Precision = precision, Recall = recall))
}

rf_pr <- get_pr(cm_rf)
svm_pr <- get_pr(cm_svm)

pr_df <- tibble(
  Model = rep(c("Random Forest", "SVM"), each = 2),
  Metric = rep(c("Precision", "Recall"), 2),
  Value = c(rf_pr, svm_pr)
)

# Bar plot
pr_plot <- ggplot(pr_df, aes(x = Metric, y = Value, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#d7191c", "#2c7bb6")) +
  ylim(0, 1) +
  theme_minimal() +
  labs(title = "Precision and Recall for Functional Class", y = "Score")
print(pr_plot)

dev.off()

# Generate confusion matrices on test data
cm_rf <- confusionMatrix(predict(rf_model, test_data), test_data$Class)
cm_svm <- confusionMatrix(predict(svm_model, test_data), test_data$Class)

# Extract precision & recall for both classes
metrics_df <- tibble(
  Class = rep(c("Functional", "NonFunctional"), each = 2),
  Metric = rep(c("Precision", "Recall"), times = 2),
  RF = c(cm_rf$byClass["Pos Pred Value"], cm_rf$byClass["Sensitivity"],
         cm_rf$byClass["Neg Pred Value"], cm_rf$byClass["Specificity"]),
  SVM = c(cm_svm$byClass["Pos Pred Value"], cm_svm$byClass["Sensitivity"],
          cm_svm$byClass["Neg Pred Value"], cm_svm$byClass["Specificity"])
)

# Reshape for plotting
plot_df <- metrics_df %>%
  pivot_longer(cols = c("RF", "SVM"), names_to = "Model", values_to = "Score")

library(ggplot2)

ggplot(plot_df, aes(x = Metric, y = Score, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ Class) +
  scale_fill_manual(values = c("red", "#2c7bb6")) +
  labs(title = "Precision and Recall by Class and Model", y = "Score", x = "Metric") +
  theme_minimal()
true_counts <- table(test_data$Class)
rf_pred_counts <- table(predict(rf_model, test_data))
svm_pred_counts <- table(predict(svm_model, test_data))

count_df <- tibble(
  Class = rep(c("Functional", "NonFunctional"), 3),
  Source = rep(c("True", "RF Prediction", "SVM Prediction"), each = 2),
  Count = c(true_counts, rf_pred_counts, svm_pred_counts)
)

ggplot(count_df, aes(x = Class, y = Count, fill = Source)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("red", "black", "#2c7bb6")) +
  labs(title = "True vs Predicted Functional Class Counts", y = "Number of Transcripts", x = "") +
  theme_minimal()
