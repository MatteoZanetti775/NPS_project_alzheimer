## =========================================================
## Standalone script (TERM 5 REMOVED):
##  1) Build CN-based deviation terms (L1, Var, Skew, Quad)
##  2) Print correlation matrix
##  3) 5-fold stratified CV classification (one term at a time)
##     -> confusion matrix + accuracy per term
##  4) Progress bars where things may take time
## =========================================================

## ---- Packages ----
pkgs <- c("dplyr", "stringr", "ggplot2", "e1071", "corpcor", "nnet", "progress")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(dplyr)
library(stringr)
library(e1071)
library(corpcor)
library(nnet)
library(progress)

set.seed(123)

## ---- 0) Load data ----
data <- read.csv("dataQR_clean.csv")
data$DIAGNOSIS <- factor(data$DIAGNOSIS, levels = c("CN","MCI","AD"))

## ---- 1) Feature columns + numeric coercion ----
feature_cols <- setdiff(colnames(data), c("RID","DIAGNOSIS","GENOTYPE","DIAGNOSIS_NUM"))
data[feature_cols] <- lapply(data[feature_cols], as.numeric)

## ---- 2) Remove extreme row outliers (optional but recommended) ----
threshold <- 1888  # keep as you used; adjust if needed

med_all <- apply(data[feature_cols], 2, median, na.rm = TRUE)
mad_all <- apply(data[feature_cols], 2, mad,    na.rm = TRUE)
mad_all[mad_all == 0] <- 1e-6

robust_z_all <- sweep(data[feature_cols], 2, med_all, "-")
robust_z_all <- sweep(robust_z_all, 2, mad_all, "/")

extreme_rows <- apply(abs(robust_z_all) > threshold, 1, any)
cat("Rows removed as extreme outliers:", sum(extreme_rows), "\n")

data_clean <- data[!extreme_rows, ]
cat("Remaining rows:", nrow(data_clean), "\n\n")

## ---- 3) Build CN-based robust Z (baseline) ----
cn_df <- data_clean %>% filter(DIAGNOSIS == "CN")

m <- apply(cn_df[, feature_cols, drop = FALSE], 2, median, na.rm = TRUE)
s <- apply(cn_df[, feature_cols, drop = FALSE], 2, mad,    na.rm = TRUE)
s[s == 0] <- 1e-6

X <- as.matrix(data_clean[, feature_cols, drop = FALSE])
Z <- sweep(X, 2, m, "-")
Z <- sweep(Z, 2, s, "/")
absZ <- abs(Z)

p <- length(feature_cols)
eps <- 1e-12

## ---- 4) Define terms (NO TERM 5) ----
## Term1: mean absolute deviation (L1), logged
term1_L1 <- rowSums(absZ, na.rm = TRUE) / p
log_term1_L1 <- log(term1_L1 + eps)

## Term3: variance of |Z|, logged
term3_var_absZ <- apply(absZ, 1, var, na.rm = TRUE)
log_term3 <- log(term3_var_absZ + eps)

## Term4_1: skewness of |Z|
term4_1_skew_absZ <- apply(absZ, 1, function(v) e1071::skewness(v, na.rm = TRUE, type = 2))

## Term4_2: correlation-pattern deviation (quadratic form)
cat("Computing Term 4_2 (quadratic form using shrinkage covariance)...\n")
pb_q <- progress_bar$new(
  format = "Term 4_2 [:bar] :percent | ETA: :eta",
  total  = 3,
  clear  = FALSE,
  width  = 60
)

Z_cn <- Z[data_clean$DIAGNOSIS == "CN", , drop = FALSE]
pb_q$tick()

Sigma_cn <- corpcor::cov.shrink(Z_cn)
pb_q$tick()

Sigma_inv <- solve(Sigma_cn)
term4_2_quad <- apply(Z, 1, function(z) as.numeric(t(z) %*% Sigma_inv %*% z))
pb_q$tick()
cat("Done.\n\n")

## ---- 5) Collect into one data.frame ----
scored <- data_clean %>%
  mutate(
    log_term1_L1         = log_term1_L1,
    log_term3            = log_term3,
    term4_1_skew_absZ    = term4_1_skew_absZ,
    term4_2_quad         = term4_2_quad
  )

## ---- 6) Correlation matrix (NO TERM 5) ----
terms_df <- scored[, c("log_term1_L1", "log_term3", "term4_1_skew_absZ", "term4_2_quad")]
cor_mat <- cor(terms_df, use = "complete.obs")
cat("Correlation matrix (complete cases):\n")
print(round(cor_mat, 4))
cat("\n")

## =========================================================
## 7) Stratified 5-fold CV (per-term, 1D multinomial logistic)
##    Output: confusion matrix + total accuracy per term
## =========================================================

K <- 5

term_names <- c("log_term1_L1", "log_term3", "term4_1_skew_absZ", "term4_2_quad")

## Progress bar: terms * folds
pb <- progress_bar$new(
  format = "CV [:bar] :current/:total (:percent) | ETA: :eta",
  total  = length(term_names) * K,
  clear  = FALSE,
  width  = 70
)

cv_results <- list()

for (tname in term_names) {
  
  # Work on complete cases for this term
  df_t <- scored %>%
    select(DIAGNOSIS, all_of(tname)) %>%
    filter(is.finite(.data[[tname]]), !is.na(DIAGNOSIS))
  
  # Rebuild fold assignment for this reduced set (keeps stratification)
  f_id <- rep(NA_integer_, nrow(df_t))
  for (lev in levels(df_t$DIAGNOSIS)) {
    idx <- which(df_t$DIAGNOSIS == lev)
    idx <- sample(idx)
    f_id[idx] <- rep(1:K, length.out = length(idx))
  }
  
  all_pred <- rep(NA_character_, nrow(df_t))
  all_true <- df_t$DIAGNOSIS
  
  for (k in 1:K) {
    test_idx  <- which(f_id == k)
    train_idx <- which(f_id != k)
    
    trainD <- df_t[train_idx, , drop = FALSE]
    testD  <- df_t[test_idx,  , drop = FALSE]
    
    # 1D multinomial logistic regression
    model <- nnet::multinom(
      as.formula(paste0("DIAGNOSIS ~ ", tname)),
      data = trainD,
      trace = FALSE
    )
    
    pred <- predict(model, newdata = testD)
    all_pred[test_idx] <- as.character(pred)
    
    pb$tick()
  }
  
  all_pred <- factor(all_pred, levels = levels(all_true))
  conf <- table(Predicted = all_pred, True = all_true)
  acc  <- mean(all_pred == all_true, na.rm = TRUE)
  
  cv_results[[tname]] <- list(accuracy = acc, confusion = conf)
  
  cat("========================================\n")
  cat("Term:", tname, "\n")
  cat("Total accuracy:", round(acc, 4), "\n")
  cat("Confusion matrix:\n")
  print(conf)
  cat("\n")
}

## ---- Optional: quick accuracy summary table ----
acc_summary <- data.frame(
  term = names(cv_results),
  accuracy = sapply(cv_results, function(x) x$accuracy)
) %>%
  arrange(desc(accuracy))

cat("========================================\n")
cat("Accuracy summary (descending):\n")
print(acc_summary)
