library(tidyverse)

#-----Read data-----
data <- read.csv("Doc_for_ANOVA.csv", header = T)

#-----Convert to factors-----
data$Condition <- factor(data$Condition)
data$gene <- factor(data$gene)

#-----perform ANOVA for all genes-----
genes <- unique(data$gene)

all_anova_results <- data.frame()
all_descriptive_stats <- data.frame()
all_tukey_results <- data.frame()
summary_results <- data.frame()

# Loop through each gene
for (gene in genes) {
  
  cat("Analyzing gene:", as.character(gene), "\n")
  
  # Subset data for this gene
  gene_data <- data %>% filter(gene == !!gene)
  
  # Check if there are at least 2 conditions
  unique_conditions <- unique(gene_data$Condition)
  if (length(unique_conditions) < 2) {
    cat("  Skipped - only one condition present\n")
    next
  }
  
  # Remove any NA or empty conditions
  gene_data <- gene_data %>% filter(!is.na(Condition) & Condition != "")
  
  # Re-level the factor to remove unused levels
  gene_data$Condition <- droplevels(gene_data$Condition)
  
  # Check again after cleaning
  if (length(levels(gene_data$Condition)) < 2) {
    cat("  Skipped - insufficient conditions after cleaning\n")
    next
  }
  
  # ---- Descriptive Statistics ----
  stats <- gene_data %>%
    group_by(Condition) %>%
    summarise(
      Mean = mean(Value),
      SD = sd(Value),
      SEM = sd(Value) / sqrt(n()),
      n = n(),
      .groups = "drop"
    )
  stats$Gene <- as.character(gene)
  all_descriptive_stats <- rbind(all_descriptive_stats, stats)
  
  # ---- ANOVA ----
  anova_model <- tryCatch({
    aov(Value ~ Condition, data = gene_data)
  }, error = function(e) {
    cat("  Error in ANOVA:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(anova_model)) {
    next
  }
  
  anova_summary <- summary(anova_model)
  anova_table <- anova_summary[[1]]
  
  # Extract ANOVA results
  p_value <- anova_table[1, 5]
  f_stat <- anova_table[1, 4]
  
  if (is.na(p_value) || is.na(f_stat)) {
    cat("  Skipped - ANOVA could not be calculated\n")
    next
  }
  
  # Store ANOVA table
  anova_df <- data.frame(
    Gene = as.character(gene),
    Source = c("Between_Groups", "Within_Groups"),
    Df = anova_table[, 1],
    Sum_Sq = anova_table[, 2],
    Mean_Sq = anova_table[, 3],
    F_value = c(f_stat, NA),
    P_value = c(p_value, NA)
  )
  all_anova_results <- rbind(all_anova_results, anova_df)
  
  # ---- Tukey's HSD ----
  tukey_result <- tryCatch({
    TukeyHSD(anova_model, conf.level = 0.95)
  }, error = function(e) {
    return(NULL)
  })
  
  if (!is.null(tukey_result)) {
    tukey_df <- as.data.frame(tukey_result$Condition)
    tukey_df$Gene <- as.character(gene)
    tukey_df$Comparison <- rownames(tukey_df)
    rownames(tukey_df) <- NULL
    all_tukey_results <- rbind(all_tukey_results, tukey_df)
  }
  
  # ---- Summary ----
  # Determine significance
  if (p_value < 0.001) {
    sig <- "***"
  } else if (p_value < 0.01) {
    sig <- "**"
  } else if (p_value < 0.05) {
    sig <- "*"
  } else {
    sig <- "ns"
  }
  
  # Create summary row
  summary_row <- data.frame(Gene = as.character(gene))
  
  # Add means and SDs for each condition
  for (cond in unique(gene_data$Condition)) {
    cond_stats <- stats[stats$Condition == cond, ]
    summary_row[[paste0("Mean_", cond)]] <- cond_stats$Mean
    summary_row[[paste0("SD_", cond)]] <- cond_stats$SD
  }
  
  summary_row$F_statistic <- f_stat
  summary_row$P_value <- p_value
  summary_row$Significance <- sig
  
  summary_results <- rbind(summary_results, summary_row)
  
  cat("  F =", round(f_stat, 2), ", P =", format(p_value, scientific = FALSE, digits = 4), sig, "\n")
}

# ============================================
# DISPLAY RESULTS
# ============================================

cat("\n==================================================\n")
cat("ANOVA SUMMARY\n")
cat("==================================================\n")
print(summary_results, row.names = FALSE)

cat("\n==================================================\n")
cat("DESCRIPTIVE STATISTICS\n")
cat("==================================================\n")
print(all_descriptive_stats, row.names = FALSE)

# ============================================
# SAVE RESULTS
# ============================================

write.csv(all_anova_results, "ANOVA_Complete_Table.csv", row.names = FALSE)
write.csv(all_descriptive_stats, "ANOVA_Descriptive_Stats.csv", row.names = FALSE)
write.csv(all_tukey_results, "ANOVA_Tukey_Results.csv", row.names = FALSE)
write.csv(summary_results, "ANOVA_Summary.csv", row.names = FALSE)

cat("\n==================================================\n")
cat("ANOVA Analysis Complete!\n")
cat("==================================================\n")
cat("Files saved:\n")
cat("  1. ANOVA_Complete_Table.csv - Full ANOVA tables\n")
cat("  2. ANOVA_Descriptive_Stats.csv - Means, SD, SEM\n")
cat("  3. ANOVA_Tukey_Results.csv - Tukey's HSD results\n")
cat("  4. ANOVA_Summary.csv - Summary with significance\n")
cat("==================================================\n")
