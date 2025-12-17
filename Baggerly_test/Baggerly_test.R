# ==============================
# SAFE BAGGERLY'S TEST FOR RNA-SEQ
# ==============================
getwd()
# 1) Load counts
counts <- read.csv("Combined_counts.csv", row.names = 1)

# Ensure all counts are numeric
counts[] <- lapply(counts, as.numeric)

# 2) Define groups
# B vs D
group <- c("Glucose", "Glucose", "Glucose", "Cellobiose", "Cellobiose", "Cellobiose")
if(length(group) != ncol(counts)) stop("Group length does not match number of columns")

# 3) Add pseudo-count to avoid zeros
counts <- counts + 1

# 4) Filter low-count genes
keep_genes <- rowSums(counts[, group == "Glucose"]) > 1 & rowSums(counts[, group == "Cellobiose"]) > 1
counts <- counts[keep_genes, ]

# 5) Compute library sizes and proportions
libsize <- colSums(counts)
prop_mat <- sweep(counts, 2, libsize, FUN="/")

# 6) Baggerly's test function
baggerly_test <- function(p1, p2, w1, w2) {
  pbar <- (w1*p1 + w2*p2) / (w1 + w2)
  var1 <- pbar * (1 - pbar)
  denom <- sqrt((var1/w1) + (var1/w2))
  if (denom == 0) return(NA)
  t_stat <- (p1 - p2) / denom
  pval <- 2 * (1 - pnorm(abs(t_stat)))
  return(pval)
}

# 7) Run test for all genes
results <- data.frame(
  gene = rownames(prop_mat),
  p.value = NA_real_,
  log2FC = NA_real_
)

for (i in 1:nrow(prop_mat)) {
  pp <- as.numeric(prop_mat[i, ])
  p1 <- mean(pp[group == "Glucose"], na.rm = TRUE)
  p2 <- mean(pp[group == "Cellobiose"], na.rm = TRUE)
  w1 <- sum(libsize[group == "Glucose"])
  w2 <- sum(libsize[group == "Cellobiose"])
  
  # log2 fold change
  results$log2FC[i] <- log2((p1 + 1e-8) / (p2 + 1e-8))
  
  # skip if pooled proportion = 0 or 1
  pbar <- (w1*p1 + w2*p2)/(w1 + w2)
  if (pbar == 0 | pbar == 1) {
    results$p.value[i] <- NA
  } else {
    results$p.value[i] <- baggerly_test(p1, p2, w1, w2)
  }
}

# 8) Adjust for multiple testing
results$adj.p.value <- p.adjust(results$p.value, method = "BH")

# 9) Sort by adjusted p-value
results <- results[order(results$adj.p.value), ]

# 10) Save results
write.csv(results, "Baggerly_test_results.csv", row.names = FALSE)

# ==============================
# END
# ==============================
