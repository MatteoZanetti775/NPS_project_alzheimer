## =========================================================
## TERM 1 ONLY (optimized loop)
## - Loads dataQR_clean.csv
## - Optional extreme-row removal (robust z vs global med/MAD)
## - Computes CN-based robust Z per feature (median/MAD on CN)
## - Term 1 = mean_j |Z_ij|  (L1 across regions), then log
## - Outputs: data frame with RID, DIAGNOSIS, GENOTYPE, log_term1_L1
## =========================================================

pkgs <- c("dplyr")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)
library(dplyr)

set.seed(123)

## ---- 0) Load data ----
data <- read.csv("dataQR_clean.csv")
data$DIAGNOSIS <- factor(data$DIAGNOSIS, levels = c("CN","MCI","AD"))

## ---- 1) Feature columns + numeric coercion ----
feature_cols <- setdiff(colnames(data), c("RID","DIAGNOSIS","GENOTYPE","DIAGNOSIS_NUM"))
data[feature_cols] <- lapply(data[feature_cols], as.numeric)

## ---- 2) Remove extreme row outliers (optional but consistent with your script) ----
threshold <- 1888  # keep same default
med_all <- apply(data[feature_cols], 2, median, na.rm = TRUE)
mad_all <- apply(data[feature_cols], 2, mad,    na.rm = TRUE)
mad_all[mad_all == 0] <- 1e-6

# robust z (global) without building huge intermediates too much
# (this part is still matrix-based; you can skip it by setting threshold <- Inf)
X_all <- as.matrix(data[feature_cols])
robust_z_all <- sweep(X_all, 2, med_all, "-")
robust_z_all <- sweep(robust_z_all, 2, mad_all, "/")

extreme_rows <- apply(abs(robust_z_all) > threshold, 1, any)
cat("Rows removed as extreme outliers:", sum(extreme_rows), "\n")

data_clean <- data[!extreme_rows, ]
cat("Remaining rows:", nrow(data_clean), "\n\n")

rm(X_all, robust_z_all)  # free memory

## ---- 3) CN-based med/MAD per feature ----
cn_df <- data_clean %>% filter(DIAGNOSIS == "CN")
m_cn <- apply(cn_df[, feature_cols, drop = FALSE], 2, median, na.rm = TRUE)
s_cn <- apply(cn_df[, feature_cols, drop = FALSE], 2, mad,    na.rm = TRUE)
s_cn[s_cn == 0 | is.na(s_cn)] <- 1e-6

## ---- 4) TERM 1 (optimized loop) ----
# Term 1 = mean_j | (x_ij - m_j) / s_j |
# Do NOT build Z or absZ; accumulate row sums in a loop.

X <- as.matrix(data_clean[, feature_cols, drop = FALSE])
n <- nrow(X)
p <- ncol(X)

abs_sum <- numeric(n)

pb <- txtProgressBar(min = 0, max = p, style = 3)
for (j in seq_len(p)) {
  zj <- (X[, j] - m_cn[j]) / s_cn[j]
  abs_sum <- abs_sum + abs(zj)
  setTxtProgressBar(pb, j)
}
close(pb)

term1_L1 <- abs_sum / p
eps <- 1e-12
log_term1_L1 <- log(term1_L1 + eps)

## ---- 5) Output (only Term 1) ----
term1_df <- data_clean %>%
  transmute(
    RID,
    DIAGNOSIS,
    GENOTYPE,
    log_term1_L1
  )

print(head(term1_df))
cat("\nComputed log_term1_L1 for", nrow(term1_df), "subjects.\n")

