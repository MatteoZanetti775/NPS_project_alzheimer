# ============================================================
# SECTION 5 (FULL SCRIPT - FIXED)
# Logistic regression - (CN+MCI) vs AD   [binary: AD vs non-AD]
# 5-fold CV on STABLE CN+MCI+AD (patient-level stratified splits)
#
# Predictors: functional_depth + age + GENOTYPE
# Depth is computed as mean MBD vs reference curves (TRAIN-only per fold)
#   - Reference = CN if present, else (CN+MCI) within TRAIN fold
#
# Outputs:
#  - Confusion matrix (CV aggregated on stable)
#  - CV accuracy (stable)
#  - CV balanced accuracy (stable)
#  - Accuracy on converters (final diagnosis: AD vs non-AD), trained on ALL stable
#
# Progress:
#  - VERY granular % progress for curve-building per fold (train + test)
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
all_volume_names <- colnames(datatot)[17:133]

# Stable CN/MCI/AD only
data_stable <- datatot %>%
  filter(nav != 1) %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) == 1) %>%
  ungroup() %>%
  mutate(DIAG = case_when(
    DIAGNOSIS == 1 ~ "CN",
    DIAGNOSIS == 2 ~ "MCI",
    DIAGNOSIS == 4 ~ "AD",
    TRUE ~ "OTHER"
  )) %>%
  filter(DIAG %in% c("CN", "MCI", "AD"))

# Converters: final label CN/MCI/AD, but outcome is AD vs nonAD
data_conv <- datatot %>%
  filter(nav != 1) %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) != 1) %>%
  ungroup() %>%
  mutate(DIAG_FINAL = case_when(
    last(DIAGNOSIS) == 1 ~ "CN",
    last(DIAGNOSIS) == 2 ~ "MCI",
    last(DIAGNOSIS) == 4 ~ "AD",
    TRUE ~ "OTHER"
  )) %>%
  filter(DIAG_FINAL %in% c("CN", "MCI", "AD"))

# Shared age grid
grid <- seq(min(datatot$age, na.rm = TRUE),
            max(datatot$age, na.rm = TRUE),
            length.out = GRID_L)

# ----------------------------
# 1) Patient-level labels + folds (stratified by stable class)
# ----------------------------
stable_patient_labels <- data_stable %>% distinct(RID, DIAG)

make_folds <- function(df_rid_lab, k = 5) {
  folds <- vector("list", k)
  for (i in 1:k) folds[[i]] <- character(0)
  
  for (cl in sort(unique(df_rid_lab$DIAG))) {
    rids <- df_rid_lab$RID[df_rid_lab$DIAG == cl]
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
  if (mf_ref$N == 0) stop("Reference mfData is empty: cannot compute depths.")
  depth_mat <- matrix(0, nrow = mf_test$N, ncol = mf_test$L)
  for (j in 1:mf_test$L) {
    depth_mat[, j] <- mbd_custom(mf_test$fDList[[j]]$values, mf_ref$fDList[[j]]$values)
  }
  rowMeans(depth_mat)
}

balanced_acc_binary <- function(true_lab, pred_lab) {
  lvls <- sort(unique(true_lab))
  recalls <- sapply(lvls, function(cl) {
    denom <- sum(true_lab == cl)
    if (denom == 0) return(NA_real_)
    sum(pred_lab == cl & true_lab == cl) / denom
  })
  mean(recalls, na.rm = TRUE)
}

build_patient_meta <- function(df_visits, patient_ids) {
  df_visits %>%
    filter(RID %in% patient_ids) %>%
    group_by(RID) %>%
    summarise(
      age = min(age, na.rm = TRUE),
      GENOTYPE = first(GENOTYPE),
      .groups = "drop"
    )
}

# ============================================================
# 4) Cross-validation on STABLE, outcome = AD vs nonAD
# ============================================================
all_true <- character(0)
all_pred <- character(0)

cat("Running Logistic regression (AD vs non-AD) with CV on stable patients...\n")
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
  
  train_pat <- unique(train_df$RID)
  test_pat  <- unique(test_df$RID)
  
  if (length(train_pat) == 0 || length(test_pat) == 0) {
    cat("WARNING: empty train or test patients in this fold. Skipping.\n")
    next
  }
  
  # Train-only scaling
  mu  <- colMeans(train_df[, all_volume_names], na.rm = TRUE)
  sdv <- apply(train_df[, all_volume_names], 2, sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- 1
  
  train_sc <- train_df
  test_sc  <- test_df
  train_sc[, all_volume_names] <- scale(train_sc[, all_volume_names], center = mu, scale = sdv)
  test_sc[,  all_volume_names] <- scale(test_sc[,  all_volume_names], center = mu, scale = sdv)
  
  # Reference patients in TRAIN: prefer CN; else (CN+MCI)
  ref_train_pat <- train_df %>% filter(DIAG == "CN") %>% distinct(RID) %>% pull(RID)
  if (length(ref_train_pat) == 0) {
    ref_train_pat <- train_df %>% filter(DIAG %in% c("CN","MCI")) %>% distinct(RID) %>% pull(RID)
    cat("WARNING: No CN in training fold; using (CN+MCI) as reference for depth.\n")
  }
  if (length(ref_train_pat) == 0) {
    cat("WARNING: No non-AD reference patients in training fold. Skipping fold.\n")
    next
  }
  
  # ----------------------------
  # TRAIN mfData + TRAIN reference subset
  # ----------------------------
  cat("\nBuilding TRAIN mfData (for depth + reference)...\n"); flush.console()
  mf_train <- mfData(
    grid,
    create_mf_list_progress(train_sc, train_pat, all_volume_names, grid, label = "TRAIN")
  )
  
  mf_ref <- mf_train[which(train_pat %in% ref_train_pat)]
  if (mf_ref$N == 0) {
    cat("WARNING: reference mfData ended empty. Skipping fold.\n")
    next
  }
  
  depth_train <- calc_multi_depth(mf_train, mf_ref)
  
  train_meta <- build_patient_meta(train_sc, train_pat)
  train_meta$functional_depth <- depth_train[match(train_meta$RID, train_pat)]
  train_meta <- train_meta %>%
    left_join(train_df %>% distinct(RID, DIAG), by = "RID") %>%
    mutate(Outcome = if_else(DIAG == "AD", 1, 0))
  
  # ----------------------------
  # TEST mfData, depth vs TRAIN reference
  # ----------------------------
  cat("\nBuilding TEST mfData (depth vs TRAIN reference)...\n"); flush.console()
  mf_test <- mfData(
    grid,
    create_mf_list_progress(test_sc, test_pat, all_volume_names, grid, label = "TEST ")
  )
  
  depth_test <- calc_multi_depth(mf_test, mf_ref)
  
  test_meta <- build_patient_meta(test_sc, test_pat)
  test_meta$functional_depth <- depth_test[match(test_meta$RID, test_pat)]
  test_meta <- test_meta %>%
    left_join(test_df %>% distinct(RID, DIAG), by = "RID") %>%
    mutate(Outcome = if_else(DIAG == "AD", 1, 0))
  
  # Align GENOTYPE levels (no empty-table surprises)
  train_meta$GENOTYPE <- factor(train_meta$GENOTYPE)
  test_meta$GENOTYPE  <- factor(test_meta$GENOTYPE, levels = levels(train_meta$GENOTYPE))
  
  test_meta <- test_meta %>%
    filter(!is.na(GENOTYPE), !is.na(functional_depth), !is.na(age), !is.na(Outcome))
  
  if (nrow(test_meta) == 0) {
    cat("WARNING: test_meta empty after alignment/cleanup. Skipping fold.\n")
    next
  }
  
  # Fit logistic regression
  model <- glm(Outcome ~ functional_depth + age + GENOTYPE,
               data = train_meta,
               family = binomial)
  
  prob <- predict(model, newdata = test_meta, type = "response")
  
  pred <- ifelse(prob > 0.5, "AD", "nonAD")
  true <- ifelse(test_meta$Outcome == 1, "AD", "nonAD")
  
  all_true <- c(all_true, as.character(true))
  all_pred <- c(all_pred, as.character(pred))
  
  cat(sprintf("\nFold %d appended %d predictions (end: %s)\n",
              f, length(true), format(Sys.time(), "%H:%M:%S")))
  flush.console()
}

cat(sprintf("\nTotal CV time: %.1f minutes\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- CV outputs (stable) ----
cat("\n==================== RESULTS: STABLE (CV aggregated) ====================\n")

if (length(all_true) == 0 || length(all_pred) == 0) {
  cat("ERROR: No predictions were collected across folds.\n")
  cat("Likely folds were skipped due to missing reference or empty test_meta.\n")
} else {
  cm_stable <- table(Actual = all_true, Predicted = all_pred)
  print(cm_stable)
  
  cv_acc  <- sum(diag(cm_stable)) / sum(cm_stable)
  cv_bacc <- balanced_acc_binary(all_true, all_pred)
  
  cat(sprintf("\nCV Accuracy (stable AD vs non-AD): %.4f\n", cv_acc))
  cat(sprintf("CV Balanced Accuracy (stable AD vs non-AD): %.4f\n", cv_bacc))
}

# ============================================================
# 5) Converters accuracy (train on ALL stable, test on converters)
# Outcome: FINAL diagnosis AD vs non-AD (CN or MCI)
# ============================================================
cat("\n==================== RESULTS: CONVERTERS ====================\n")

stable_all <- data_stable
conv_all   <- data_conv

stable_pat_all <- unique(stable_all$RID)
conv_pat_all   <- unique(conv_all$RID)

# Scaling from ALL stable
mu_all <- colMeans(stable_all[, all_volume_names], na.rm = TRUE)
sd_all <- apply(stable_all[, all_volume_names], 2, sd, na.rm = TRUE)
sd_all[sd_all == 0 | is.na(sd_all)] <- 1

stable_sc <- stable_all
conv_sc   <- conv_all
stable_sc[, all_volume_names] <- scale(stable_sc[, all_volume_names], center = mu_all, scale = sd_all)
conv_sc[,   all_volume_names] <- scale(conv_sc[,   all_volume_names], center = mu_all, scale = sd_all)

# Reference in ALL stable: prefer CN else (CN+MCI)
ref_all_pat <- stable_all %>% filter(DIAG == "CN") %>% distinct(RID) %>% pull(RID)
if (length(ref_all_pat) == 0) {
  ref_all_pat <- stable_all %>% filter(DIAG %in% c("CN","MCI")) %>% distinct(RID) %>% pull(RID)
  cat("WARNING: No CN in all-stable set; using (CN+MCI) as reference for depth.\n")
}

if (length(ref_all_pat) == 0) {
  cat("ERROR: No reference (non-AD) patients in all stable; cannot evaluate converters.\n")
} else {
  
  cat("Building ALL-STABLE mfData (for depth + reference)...\n"); flush.console()
  mf_stable_all <- mfData(
    grid,
    create_mf_list_progress(stable_sc, stable_pat_all, all_volume_names, grid, label = "ALLST")
  )
  
  mf_ref_all <- mf_stable_all[which(stable_pat_all %in% ref_all_pat)]
  if (mf_ref_all$N == 0) {
    cat("ERROR: reference mfData for all-stable ended empty.\n")
  } else {
    
    depth_stable_all <- calc_multi_depth(mf_stable_all, mf_ref_all)
    
    stable_meta_all <- build_patient_meta(stable_sc, stable_pat_all)
    stable_meta_all$functional_depth <- depth_stable_all[match(stable_meta_all$RID, stable_pat_all)]
    stable_meta_all <- stable_meta_all %>%
      left_join(stable_all %>% distinct(RID, DIAG), by = "RID") %>%
      mutate(Outcome = if_else(DIAG == "AD", 1, 0)) %>%
      filter(!is.na(functional_depth), !is.na(age), !is.na(GENOTYPE), !is.na(Outcome))
    
    cat("Building CONVERTERS mfData (depth vs ALL-stable reference)...\n"); flush.console()
    mf_conv <- mfData(
      grid,
      create_mf_list_progress(conv_sc, conv_pat_all, all_volume_names, grid, label = "CONV ")
    )
    
    depth_conv <- calc_multi_depth(mf_conv, mf_ref_all)
    
    conv_meta <- build_patient_meta(conv_sc, conv_pat_all)
    conv_meta$functional_depth <- depth_conv[match(conv_meta$RID, conv_pat_all)]
    conv_meta <- conv_meta %>%
      left_join(conv_all %>% distinct(RID, DIAG_FINAL), by = "RID") %>%
      mutate(Outcome = if_else(DIAG_FINAL == "AD", 1, 0)) %>%
      filter(!is.na(functional_depth), !is.na(age), !is.na(GENOTYPE), !is.na(Outcome))
    
    # Align GENOTYPE levels
    stable_meta_all$GENOTYPE <- factor(stable_meta_all$GENOTYPE)
    conv_meta$GENOTYPE <- factor(conv_meta$GENOTYPE, levels = levels(stable_meta_all$GENOTYPE))
    conv_meta <- conv_meta %>% filter(!is.na(GENOTYPE))
    
    if (nrow(conv_meta) == 0) {
      cat("ERROR: conv_meta empty after GENOTYPE alignment.\n")
    } else {
      
      model_all <- glm(Outcome ~ functional_depth + age + GENOTYPE,
                       data = stable_meta_all,
                       family = binomial)
      
      prob_c <- predict(model_all, newdata = conv_meta, type = "response")
      pred_c <- ifelse(prob_c > 0.5, "AD", "nonAD")
      true_c <- ifelse(conv_meta$Outcome == 1, "AD", "nonAD")
      
      acc_conv <- mean(pred_c == true_c)
      cat(sprintf("Accuracy on converters (final diagnosis AD vs non-AD): %.4f\n", acc_conv))
    }
  }
}

cat("=======================================================================\n")

