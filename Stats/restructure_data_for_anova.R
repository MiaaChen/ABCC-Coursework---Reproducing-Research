getwd()
library(tidyr)
library(dplyr)
library(broom)
library(multcomp)
library(ggplot2)

data <- read.csv("./Combined_RPKM.csv",header=T, stringsAsFactors = F)
colnames(data)
df_test <- data[data$gene == "AAC1", , drop = FALSE]
test_long <- df_test %>%
  pivot_longer(
    cols = -gene,
    names_to = "Sample",
    values_to = "Value"
  )

# Reshape all genes 
long_data <- data %>%
  pivot_longer(
  cols = -gene,
  names_to = "Sample",
  values_to = "Value"
)

# Parse metadata from sample names
long_data <- long_data %>%
  mutate(
    Strain    = substr(Sample, 1, 4),  # JCYL / JCTL
    Replicate = substr(Sample, 5, 7),  # 001 / 002 / 003
    Condition = substr(Sample, 8, 8)   # B / D
  )

# Convert to factors
long_data <- long_data %>%
  mutate(
    gene = factor(gene),
    Strain = factor(Strain),
    Replicate = factor(Replicate),
    Condition = factor(Condition)
  )

write.csv(long_data, "Doc_for_ANOVA.csv", row.names = F)



