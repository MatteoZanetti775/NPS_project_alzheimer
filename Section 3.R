# ============================================================
# SECTION 3 (FULL SCRIPT)
# Max-depth classifier - 2-class (CN vs AD) - ALL 117 volumes
# 5-fold CV on STABLE CN+AD only (patient-level stratified splits)
#
# Outputs:
#  - Confusion matrix (CV aggregated on stable CN/AD)
#  - CV accuracy (stable CN/AD)
#  - CV balanced accuracy (stable CN/AD)
#  - Accuracy on converters (final diagnosis CN vs AD), trained on ALL stable CN/AD
#
# Progress:
#  - VERY granular % progress for curve-building (train + test) inside each fold
# ============================================================

rm(list = ls())

library(dplyr)
library(tidyr)
library(roahd)

set.seed(123)

K_FOLDS <- 5
GRID_L  <- 50

# ----------------------------
# 0) Load + prepare data
# ----------------------------
datatot <- read.csv("data.csv")

# all 117 volumes (columns 17:133)
volume_names <- colnames(datatot)[17:133]

# Stable patients only: diagnosis never changes across visits
# Keep ONLY CN and AD here
data_stable <- datatot %>%
  filter(nav != 1) %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) == 1) %>%
  ungroup() %>%
  mutate(DIAGNOSIS = case_when(
    DIAGNOSIS == 1 ~ "CN",
    DIAGNOSIS == 4 ~ "AD",
    TRUE ~ "OTHER"
  )) %>%
  filter(DIAGNOSIS %in% c("CN", "AD"))

# Converters: diagnosis changes over time, use FINAL diagnosis as target
# Keep ONLY CN vs AD final label (drop others)
data_conv <- datatot %>%
  filter(nav != 1) %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) != 1) %>%
  ungroup() %>%
  mutate(DIAG_FINAL = case_when(
    last(DIAGNOSIS) == 1 ~ "CN",
    last(DIAGNOSIS) == 4 ~ "AD",
    TRUE ~ "OTHER"
  )) %>%
  filter(DIAG_FINAL %in% c("CN", "AD"))

# Age grid shared everywhere
grid <- seq(min(datatot$age, na.rm = TRUE),
            max(datatot$age, na.rm = TRUE),
            length.out = GRID_L)

# ----------------------------
# 1) Stratified folds by RID
# ----------------------------
stable_patient_labels <- data_stable %>% distinct(RID, DIAGNOSIS)

make_folds <- function(df_rid_lab, k = 5) {
  folds <- vector("list", k)
  for (i in 1:k) folds[[i]] <- character(0)
  
  for (cl in sort(unique(df_rid_lab$DIAGNOSIS))) {
    rids <- df_rid_lab$RID[df_rid_lab$DIAGNOSIS == cl]
    rids <- sample(rids)
    parts <- split(rids, rep(1:k, length.out = length(rids)))
    for (i in 1:k) folds[[i]] <- c(folds[[i]], parts[[i]])
  }
  folds
}

folds <- make_folds(stable_patient_labels, K_FOLDS)

# ----------------------------
# 2) VERY granular progress curve builder
# ----------------------------
create_mf_list_progress <- function(df, patient_list, volume_names, grid_points, label = "") {
  
  total_steps <- length(patient_list) * length(volume_names)
  step_counter <- 0
  t_start <- Sys.time()
  
  # Pre-split by RID to avoid repeated filtering
  df_by_patient <- split(df, df$RID)
  
  update_progress <- function(step) {
    pct <- 100 * step / total_steps
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    rate <- step / max(elapsed, 1e-6)
    remaining <- (total_steps - step) / max(rate, 1e-6)
    
    cat(sprintf(
      "\r%s Progress: %6.2f%% | Step %d/%d | Elapsed: %.1fs | ETA: %.1fs",
      label, pct, step, total_steps, elapsed, remaining
    ))
    flush.console()
  }
  
  result <- lapply(volume_names, function(vol) {
    
    mat <- t(sapply(patient_list, function(p) {
      
      p_data <- df_by_patient[[as.character(p)]]
      p_ag <- aggregate(p_data[[vol]] ~ p_data$age, FUN = mean)
      
      y <- approx(p_ag[[1]], p_ag[[2]], xout = grid_points, rule = 2)$y
      
      step_counter <<- step_counter + 1
      update_progress(step_counter)
      
      y
    }))
    
    mat
  })
  
  cat("\n")
  return(result)
}

# ----------------------------
# 3) Custom MBD + multi-depth
# ----------------------------
mbd_custom <- function(X, P) {
  if (is.vector(X)) X <- matrix(X, nrow = 1)
  n <- nrow(P); L <- ncol(P)
  depths <- numeric(nrow(X))
  for (i in 1:nrow(X)) {
    count <- 0
    for (t in 1:L) {
      n_a <- sum(P[, t] <= X[i, t])
      n_b <- sum(P[, t] >= X[i, t])
      count <- count + (n_a * n_b)
    }
    depths[i] <- count / (L * n * (n + 1) / 2)
  }
  depths
}

calc_multi_depth <- function(mf_test, mf_ref) {
  depth_mat <- matrix(0, nrow = mf_test$N, ncol = mf_test$L)
  for (j in 1:mf_test$L) {
    depth_mat[, j] <- mbd_custom(mf_test$fDList[[j]]$values, mf_ref$fDList[[j]]$values)
  }
  rowMeans(depth_mat)
}

balanced_acc_binary <- function(true_lab, pred_lab) {
  # mean(recall_CN, recall_AD)
  lvls <- sort(unique(true_lab))
  recalls <- sapply(lvls, function(cl) {
    denom <- sum(true_lab == cl)
    if (denom == 0) return(NA_real_)
    sum(pred_lab == cl & true_lab == cl) / denom
  })
  mean(recalls, na.rm = TRUE)
}

# ============================================================
# 4) Cross-validation on STABLE CN/AD
# ============================================================
all_true <- character(0)
all_pred <- character(0)

cat("Running 2-class (CN vs AD) max-depth CV on stable patients (ALL volumes)...\n")
t0 <- Sys.time()

for (f in 1:K_FOLDS) {
  
  cat("\n============================================================\n")
  cat(sprintf("FOLD %d/%d (start: %s)\n", f, K_FOLDS, format(Sys.time(), "%H:%M:%S")))
  cat("============================================================\n")
  flush.console()
  
  test_rids  <- folds[[f]]
  train_rids <- setdiff(stable_patient_labels$RID, test_rids)
  
  train_df <- data_stable %>% filter(RID %in% train_rids)
  test_df  <- data_stable %>% filter(RID %in% test_rids)
  
  # Scale using TRAIN only (no leakage)
  mu  <- colMeans(train_df[, volume_names], na.rm = TRUE)
  sdv <- apply(train_df[, volume_names], 2, sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- 1
  
  train_df[, volume_names] <- scale(train_df[, volume_names], center = mu, scale = sdv)
  test_df[,  volume_names] <- scale(test_df[,  volume_names], center = mu, scale = sdv)
  
  train_pat <- unique(train_df$RID)
  test_pat  <- unique(test_df$RID)
  
  cat("\nBuilding TRAIN curves (mfData)...\n"); flush.console()
  mf_train <- mfData(grid, create_mf_list_progress(train_df, train_pat, volume_names, grid, label = "TRAIN"))
  
  cat("\nBuilding TEST curves (mfData)...\n"); flush.console()
  mf_test  <- mfData(grid, create_mf_list_progress(test_df,  test_pat,  volume_names, grid, label = "TEST "))
  
  # Class references from TRAIN only
  cn_rids <- train_df %>% filter(DIAGNOSIS == "CN") %>% distinct(RID) %>% pull(RID)
  ad_rids <- train_df %>% filter(DIAGNOSIS == "AD") %>% distinct(RID) %>% pull(RID)
  
  mf_cn <- mf_train[which(train_pat %in% cn_rids)]
  mf_ad <- mf_train[which(train_pat %in% ad_rids)]
  
  cat("\nComputing depths for TEST patients...\n"); flush.console()
  d_cn <- calc_multi_depth(mf_test, mf_cn)
  d_ad <- calc_multi_depth(mf_test, mf_ad)
  
  pred <- ifelse(d_cn >= d_ad, "CN", "AD")
  
  true <- test_df %>%
    distinct(RID, DIAGNOSIS) %>%
    arrange(match(RID, test_pat)) %>%
    pull(DIAGNOSIS)
  
  all_true <- c(all_true, true)
  all_pred <- c(all_pred, pred)
  
  cat(sprintf("\nFold %d done (end: %s)\n", f, format(Sys.time(), "%H:%M:%S")))
  flush.console()
}

cat(sprintf("\nTotal CV time: %.1f minutes\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- CV outputs (stable) ----
cat("\n==================== RESULTS: STABLE (CV aggregated) ====================\n")
cm_stable <- table(Actual = all_true, Predicted = all_pred)
print(cm_stable)

cv_acc  <- sum(diag(cm_stable)) / sum(cm_stable)
cv_bacc <- balanced_acc_binary(all_true, all_pred)

cat(sprintf("\nCV Accuracy (stable CN/AD): %.4f\n", cv_acc))
cat(sprintf("CV Balanced Accuracy (stable CN/AD): %.4f\n", cv_bacc))

# ============================================================
# 5) Converters accuracy (train on ALL stable CN/AD, test on converters CN/AD)
# ============================================================
cat("\n==================== RESULTS: CONVERTERS ====================\n")

# scale converters using ALL stable stats (stable CN/AD only)
mu_all <- colMeans(data_stable[, volume_names], na.rm = TRUE)
sd_all <- apply(data_stable[, volume_names], 2, sd, na.rm = TRUE)
sd_all[sd_all == 0 | is.na(sd_all)] <- 1

stable_sc <- data_stable
conv_sc   <- data_conv

stable_sc[, volume_names] <- scale(stable_sc[, volume_names], center = mu_all, scale = sd_all)
conv_sc[,   volume_names] <- scale(conv_sc[,   volume_names], center = mu_all, scale = sd_all)

stable_pat_all <- unique(stable_sc$RID)
conv_pat_all   <- unique(conv_sc$RID)

cat("Building ALL stable curves...\n"); flush.console()
mf_stable_all <- mfData(grid, create_mf_list_progress(stable_sc, stable_pat_all, volume_names, grid, label = "ALLST"))

cat("Building converters curves...\n"); flush.console()
mf_conv_all   <- mfData(grid, create_mf_list_progress(conv_sc, conv_pat_all, volume_names, grid, label = "CONV "))

# references from ALL stable
cn_rids_all <- stable_sc %>% filter(DIAGNOSIS == "CN") %>% distinct(RID) %>% pull(RID)
ad_rids_all <- stable_sc %>% filter(DIAGNOSIS == "AD") %>% distinct(RID) %>% pull(RID)

mf_cn_all <- mf_stable_all[which(stable_pat_all %in% cn_rids_all)]
mf_ad_all <- mf_stable_all[which(stable_pat_all %in% ad_rids_all)]

cat("Predicting converters (final diagnosis CN vs AD)...\n"); flush.console()
d_cn_c <- calc_multi_depth(mf_conv_all, mf_cn_all)
d_ad_c <- calc_multi_depth(mf_conv_all, mf_ad_all)

pred_conv <- ifelse(d_cn_c >= d_ad_c, "CN", "AD")

true_conv <- conv_sc %>%
  group_by(RID) %>%
  summarise(DIAG_FINAL = first(DIAG_FINAL), .groups = "drop") %>%
  arrange(match(RID, conv_pat_all)) %>%
  pull(DIAG_FINAL)

acc_conv <- mean(pred_conv == true_conv)

cat(sprintf("Accuracy on converters (final diagnosis CN vs AD): %.4f\n", acc_conv))
cat("=======================================================================\n")
