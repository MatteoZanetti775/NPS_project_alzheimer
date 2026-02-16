# ============================================================
# Functional CN-envelope ordinal deviation classifier (Term 1)
# 3-class: CN vs MCI vs AD
# Patient-level stratified 5-fold CV
#
# Core idea:
#  - Build CN median/MAD envelopes over time (per region, per time grid point)
#  - For each patient, count how many (region,time) points fall:
#       CN band -> 0, Ring -> w1, Outside -> w2
#  - Score = (w1*n_ring + w2*n_out) / n_total
#  - Tune (w1,w2) inside each fold (FAST: precompute counts once)
#  - Tune 2 thresholds (t1,t2) on TRAIN for max accuracy
#
# Progress bars:
#  - Envelope building (per region)
#  - Tier-count precompute (per patient)
#  - (w1,w2) search (coarse + refine)
#  - Threshold search (t1,t2)
#  - Test scoring (per patient)
#
# Notes:
#  - Uses <=4 earliest visits per patient
#  - Uses percent-change normalization per patient by default
#  - Optionally subsamples regions for tuning to reduce load
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(progress)
})

set.seed(123)

# ----------------------------
# USER SETTINGS
# ----------------------------
K_FOLDS        <- 5
max_visits     <- 4
region_start   <- 17
classes        <- c(1, 2, 4)      # CN, MCI, AD
grid_step      <- 6              # months
Tmax_use       <- 48             # cap horizon
use_pct_change <- TRUE           # recommended

# Envelope widths (fixed)
a <- 1.0                          # inner CN band multiplier
b <- 2.5                          # outer band multiplier (ring is between a and b)

# -------- Penalty tuning grids (FAST) --------
# Coarse search (small)
w1_vals_coarse <- seq(0, 10, by = 1.0)
w2_vals_coarse <- seq(0, 15, by = 1.0)

# Refine around best (still fast)
w1_refine_window <- 1.5
w2_refine_window <- 2.0
w1_refine_step   <- 0.25
w2_refine_step   <- 0.25

# Threshold tuning grid size (quantile grid)
n_thr_grid <- 80

# -------- Optional: reduce regions for tuning only (keeps final scoring on ALL regions) --------
# If you have tons of regions, tuning can still be slower due to tier-count precompute.
# This reduces that cost while keeping final chosen (w1,w2) applied to all regions.
use_region_subset_for_tuning <- TRUE
regions_tune_max <- 40   # e.g. 30-60 is reasonable

# ----------------------------
# 0) Load + basic cleanup
# ----------------------------
dat <- read.csv("data.csv")
dat$DIAGNOSIS <- as.numeric(as.character(dat$DIAGNOSIS))
dat$month     <- as.numeric(as.character(dat$month))

region_cols  <- region_start:ncol(dat)
region_names <- colnames(dat)[region_cols]

df0 <- dat %>%
  filter(!is.na(RID), !is.na(month), !is.na(DIAGNOSIS)) %>%
  filter(DIAGNOSIS %in% classes) %>%
  select(RID, month, DIAGNOSIS, all_of(region_names))

# Ensure unique names
colnames(df0) <- make.unique(colnames(df0))
region_names  <- colnames(df0)[4:ncol(df0)]

label_class <- function(x) factor(x, levels = c(1,2,4), labels = c("CN","MCI","AD"))

# ----------------------------
# Helper: patient-level stratified folds
# ----------------------------
patient_labels <- df0 %>% distinct(RID, DIAGNOSIS)

make_folds <- function(df_rid_lab, k = 5) {
  folds <- vector("list", k)
  for (i in 1:k) folds[[i]] <- integer(0)
  
  for (cl in sort(unique(df_rid_lab$DIAGNOSIS))) {
    rids <- df_rid_lab$RID[df_rid_lab$DIAGNOSIS == cl]
    rids <- sample(rids)
    parts <- split(rids, rep(1:k, length.out = length(rids)))
    for (i in 1:k) folds[[i]] <- c(folds[[i]], parts[[i]])
  }
  folds
}

folds <- make_folds(patient_labels, K_FOLDS)

# ----------------------------
# Helper: interpolation to grid
# ----------------------------
interp_to_grid <- function(t, y, grid) {
  if (length(t) < 2) return(rep(NA_real_, length(grid)))
  approx(x = t, y = y, xout = grid, method = "linear", rule = 1, ties = "ordered")$y
}

# ----------------------------
# Helper: build patient curves matrix for one region (n_pat x n_grid)
# using <= max_visits visits per patient
# ----------------------------
build_region_curves <- function(df, patient_ids, region, grid,
                                max_visits = 4, use_pct_change = TRUE) {
  
  df_sub <- df %>%
    filter(RID %in% patient_ids) %>%
    select(RID, month, all_of(region)) %>%
    rename(y = all_of(region)) %>%
    filter(!is.na(y)) %>%
    arrange(RID, month) %>%
    group_by(RID) %>%
    slice_head(n = max_visits) %>%
    ungroup()
  
  out <- matrix(NA_real_, nrow = length(patient_ids), ncol = length(grid))
  rownames(out) <- as.character(patient_ids)
  
  split_pat <- split(df_sub, df_sub$RID)
  
  for (i in seq_along(patient_ids)) {
    rid <- patient_ids[i]
    si <- split_pat[[as.character(rid)]]
    if (is.null(si) || nrow(si) < 2) next
    
    t <- si$month
    y <- si$y
    
    if (use_pct_change) {
      if (!is.finite(y[1]) || y[1] == 0) next
      y <- 100 * (y - y[1]) / y[1]
    }
    
    out[i, ] <- interp_to_grid(t, y, grid)
  }
  
  out
}

# ----------------------------
# Helper: build CN envelopes for one region
# Returns list(m = median_t, s = MAD_t) each length n_grid
# ----------------------------
build_cn_envelope_region <- function(curves_cn_mat) {
  m <- apply(curves_cn_mat, 2, median, na.rm = TRUE)
  s <- apply(curves_cn_mat, 2, mad,    na.rm = TRUE)
  s[!is.finite(s) | s == 0] <- 1e-6
  list(m = m, s = s)
}

# ----------------------------
# Helper: for ONE patient, compute counts of points in CN / ring / outside
# across ALL regions and time grid (using env_list)
# Returns c(n_cn, n_ring, n_out, n_tot)
# ----------------------------
count_patient_tiers <- function(df, rid, regions, grid, env_list,
                                a, b,
                                max_visits = 4, use_pct_change = TRUE) {
  
  subp <- df %>%
    filter(RID == rid) %>%
    arrange(month) %>%
    slice_head(n = max_visits)
  
  if (nrow(subp) < 2) return(c(n_cn=0L, n_ring=0L, n_out=0L, n_tot=0L))
  
  n_cn <- 0L; n_ring <- 0L; n_out <- 0L; n_tot <- 0L
  tvec_all <- subp$month
  
  for (reg in regions) {
    yvec <- subp[[reg]]
    ok <- which(!is.na(yvec))
    if (length(ok) < 2) next
    
    t <- tvec_all[ok]
    y <- yvec[ok]
    
    if (use_pct_change) {
      if (!is.finite(y[1]) || y[1] == 0) next
      y <- 100 * (y - y[1]) / y[1]
    }
    
    y_grid <- interp_to_grid(t, y, grid)
    good <- is.finite(y_grid)
    if (!any(good)) next
    
    env <- env_list[[reg]]
    m <- env$m; s <- env$s
    
    cn_low   <- m - a * s
    cn_high  <- m + a * s
    out_low  <- m - b * s
    out_high <- m + b * s
    
    inside_cn  <- good & (y_grid >= cn_low)  & (y_grid <= cn_high)
    inside_out <- good & (y_grid >= out_low) & (y_grid <= out_high)
    
    ring <- inside_out & !inside_cn
    out  <- good & !inside_out
    
    n_cn   <- n_cn   + sum(inside_cn)
    n_ring <- n_ring + sum(ring)
    n_out  <- n_out  + sum(out)
    n_tot  <- n_tot  + sum(good)
  }
  
  c(n_cn = n_cn, n_ring = n_ring, n_out = n_out, n_tot = n_tot)
}

# ----------------------------
# Helper: score from counts (vectorized-friendly)
# ----------------------------
score_from_counts <- function(n_ring, n_out, n_tot, w1, w2) {
  sc <- (w1 * n_ring + w2 * n_out) / n_tot
  sc[!is.finite(sc)] <- NA_real_
  sc
}

# ----------------------------
# Helper: thresholding for 3-class from a 1D score
# pred = CN if score <= t1; MCI if t1 < score <= t2; AD if score > t2
# ----------------------------
predict_3class_threshold <- function(score, t1, t2) {
  ifelse(score <= t1, 1,
         ifelse(score <= t2, 2, 4))
}

# ============================================================
# CV loop
# ============================================================
all_true <- integer(0)
all_pred <- integer(0)

cat("Running functional CN-envelope Term1 classifier (3-class) ...\n")
t0 <- Sys.time()

for (f in 1:K_FOLDS) {
  
  cat("\n============================================================\n")
  cat(sprintf("FOLD %d/%d (start: %s)\n", f, K_FOLDS, format(Sys.time(), "%H:%M:%S")))
  cat("============================================================\n")
  
  test_rids  <- folds[[f]]
  train_rids <- setdiff(patient_labels$RID, test_rids)
  
  train_df <- df0 %>% filter(RID %in% train_rids)
  test_df  <- df0 %>% filter(RID %in% test_rids)
  
  train_pat <- sort(unique(train_df$RID))
  test_pat  <- sort(unique(test_df$RID))
  
  # time grid from training only
  train_first4 <- train_df %>%
    arrange(RID, month) %>%
    group_by(RID) %>%
    slice_head(n = max_visits) %>%
    ungroup()
  
  Tmax <- min(Tmax_use, max(train_first4$month, na.rm = TRUE))
  grid <- seq(0, Tmax, by = grid_step)
  cat("Grid: 0 to", Tmax, "by", grid_step, "months (n =", length(grid), ")\n")
  
  # CN patient ids in training
  cn_train_pat <- train_df %>% filter(DIAGNOSIS == 1) %>% distinct(RID) %>% pull(RID)
  if (length(cn_train_pat) < 5) {
    stop("Too few CN patients in training fold to build functional envelopes robustly.")
  }
  
  # Labels per patient (TRAIN + TEST)
  y_train <- train_df %>%
    distinct(RID, DIAGNOSIS) %>%
    arrange(match(RID, train_pat)) %>%
    pull(DIAGNOSIS)
  
  y_test <- test_df %>%
    distinct(RID, DIAGNOSIS) %>%
    arrange(match(RID, test_pat)) %>%
    pull(DIAGNOSIS)
  
  # 1) Build CN envelopes per region (TRAIN only)
  cat("Building CN envelopes per region (TRAIN only)...\n")
  
  pb_env <- progress_bar$new(
    format = "Envelopes [:bar] :current/:total (:percent) | ETA: :eta",
    total  = length(region_names),
    clear  = FALSE,
    width  = 70
  )
  
  env_list <- vector("list", length(region_names))
  names(env_list) <- region_names
  
  for (reg in region_names) {
    cn_curves <- build_region_curves(train_df, cn_train_pat, reg, grid,
                                     max_visits = max_visits,
                                     use_pct_change = use_pct_change)
    env_list[[reg]] <- build_cn_envelope_region(cn_curves)
    pb_env$tick()
  }
  
  # Choose regions used for tuning (optional)
  if (use_region_subset_for_tuning && length(region_names) > regions_tune_max) {
    set.seed(1000 + f)  # fold-stable
    regions_tune <- sample(region_names, regions_tune_max)
    cat(sprintf("\nTuning penalties using a subset of %d/%d regions.\n",
                length(regions_tune), length(region_names)))
  } else {
    regions_tune <- region_names
    cat(sprintf("\nTuning penalties using ALL %d regions.\n", length(region_names)))
  }
  
  # 2) FAST tuning (w1,w2): precompute tier counts once, then cheap grid search
  cat("\nPrecomputing tier counts for TRAIN patients (once)...\n")
  pb_cnt <- progress_bar$new(
    format = "Tier counts [:bar] :current/:total (:percent) | ETA: :eta",
    total  = length(train_pat),
    clear  = FALSE,
    width  = 70
  )
  
  cnt_mat <- matrix(0, nrow = length(train_pat), ncol = 4)
  colnames(cnt_mat) <- c("n_cn","n_ring","n_out","n_tot")
  
  for (j in seq_along(train_pat)) {
    cnt_mat[j, ] <- count_patient_tiers(
      df = train_df, rid = train_pat[j],
      regions = regions_tune, grid = grid,
      env_list = env_list,
      a = a, b = b,
      max_visits = max_visits, use_pct_change = use_pct_change
    )
    pb_cnt$tick()
  }
  
  n_ring <- cnt_mat[, "n_ring"]
  n_out  <- cnt_mat[, "n_out"]
  n_tot  <- cnt_mat[, "n_tot"]
  
  ok_pat <- is.finite(y_train) & (n_tot > 0)
  if (sum(ok_pat) < 10) stop("Too few TRAIN patients with valid tier counts for tuning.")
  
  yy <- y_train[ok_pat]
  nr <- n_ring[ok_pat]
  no <- n_out[ok_pat]
  nt <- n_tot[ok_pat]
  
  # ---- Coarse grid search (w2 > w1 only) ----
  cat("\nTuning penalties (coarse grid, w2>w1)...\n")
  grid_w <- expand.grid(w1 = w1_vals_coarse, w2 = w2_vals_coarse) %>%
    filter(w2 > w1)
  
  pb_w <- progress_bar$new(
    format = "Penalty search (coarse) [:bar] :current/:total (:percent) | ETA: :eta",
    total  = nrow(grid_w),
    clear  = FALSE,
    width  = 70
  )
  
  best_sep <- -Inf
  best_w1  <- NA_real_
  best_w2  <- NA_real_
  
  for (i in seq_len(nrow(grid_w))) {
    w1 <- grid_w$w1[i]
    w2 <- grid_w$w2[i]
    
    sc <- (w1 * nr + w2 * no) / nt
    
    group_means  <- tapply(sc, yy, mean, na.rm = TRUE)
    overall_mean <- mean(sc, na.rm = TRUE)
    n_g <- table(yy)
    
    between_var <- sum(n_g * (group_means - overall_mean)^2, na.rm = TRUE)
    within_var  <- sum(tapply(sc, yy, var, na.rm = TRUE), na.rm = TRUE)
    
    sep <- between_var / (within_var + 1e-12)
    
    if (is.finite(sep) && sep > best_sep) {
      best_sep <- sep
      best_w1  <- w1
      best_w2  <- w2
    }
    pb_w$tick()
  }
  
  cat(sprintf("\nBest coarse penalties: w1=%.3f, w2=%.3f (sep=%.4f)\n",
              best_w1, best_w2, best_sep))
  
  # ---- Refine search around best ----
  cat("\nRefining penalties around best (still w2>w1)...\n")
  w1_lo <- max(0,  best_w1 - w1_refine_window)
  w1_hi <- min(max(w1_vals_coarse), best_w1 + w1_refine_window)
  w2_lo <- max(0,  best_w2 - w2_refine_window)
  w2_hi <- min(max(w2_vals_coarse), best_w2 + w2_refine_window)
  
  w1_vals_fine <- seq(w1_lo, w1_hi, by = w1_refine_step)
  w2_vals_fine <- seq(w2_lo, w2_hi, by = w2_refine_step)
  
  grid_w2 <- expand.grid(w1 = w1_vals_fine, w2 = w2_vals_fine) %>%
    filter(w2 > w1)
  
  pb_w2 <- progress_bar$new(
    format = "Penalty search (refine) [:bar] :current/:total (:percent) | ETA: :eta",
    total  = nrow(grid_w2),
    clear  = FALSE,
    width  = 70
  )
  
  best_sep2 <- -Inf
  best_w1_2 <- best_w1
  best_w2_2 <- best_w2
  
  for (i in seq_len(nrow(grid_w2))) {
    w1 <- grid_w2$w1[i]
    w2 <- grid_w2$w2[i]
    
    sc <- (w1 * nr + w2 * no) / nt
    
    group_means  <- tapply(sc, yy, mean, na.rm = TRUE)
    overall_mean <- mean(sc, na.rm = TRUE)
    n_g <- table(yy)
    
    between_var <- sum(n_g * (group_means - overall_mean)^2, na.rm = TRUE)
    within_var  <- sum(tapply(sc, yy, var, na.rm = TRUE), na.rm = TRUE)
    
    sep <- between_var / (within_var + 1e-12)
    
    if (is.finite(sep) && sep > best_sep2) {
      best_sep2 <- sep
      best_w1_2 <- w1
      best_w2_2 <- w2
    }
    pb_w2$tick()
  }
  
  best_w1 <- best_w1_2
  best_w2 <- best_w2_2
  best_sep <- best_sep2
  
  cat(sprintf("\nBest penalties in fold %d: w1=%.3f, w2=%.3f (sep=%.4f)\n",
              f, best_w1, best_w2, best_sep))
  
  # 3) Score TRAIN patients with best (w1,w2) using ALL regions (fresh counts)
  cat("\nScoring TRAIN patients with best (w1,w2) using ALL regions...\n")
  pb_trs <- progress_bar$new(
    format = "Train scoring [:bar] :current/:total (:percent) | ETA: :eta",
    total  = length(train_pat),
    clear  = FALSE,
    width  = 70
  )
  
  score_train <- rep(NA_real_, length(train_pat))
  
  for (j in seq_along(train_pat)) {
    cnt <- count_patient_tiers(
      df = train_df, rid = train_pat[j],
      regions = region_names, grid = grid,
      env_list = env_list,
      a = a, b = b,
      max_visits = max_visits, use_pct_change = use_pct_change
    )
    score_train[j] <- score_from_counts(cnt["n_ring"], cnt["n_out"], cnt["n_tot"], best_w1, best_w2)
    pb_trs$tick()
  }
  
  # 4) Tune thresholds (t1,t2) on TRAIN to maximize accuracy
  cat("\nTuning thresholds (t1,t2) on TRAIN for max accuracy...\n")
  
  ok <- is.finite(score_train)
  sc_ok <- score_train[ok]
  yy_ok <- y_train[ok]
  
  qs <- seq(0.02, 0.98, length.out = n_thr_grid)
  thr_grid <- as.numeric(quantile(sc_ok, probs = qs, na.rm = TRUE, names = FALSE))
  thr_grid <- sort(unique(thr_grid))
  
  if (length(thr_grid) < 5) stop("Too few distinct training scores to tune thresholds.")
  
  pb_thr <- progress_bar$new(
    format = "Threshold search [:bar] :current/:total (:percent) | ETA: :eta",
    total  = length(thr_grid) * length(thr_grid),
    clear  = FALSE,
    width  = 70
  )
  
  best_acc <- -Inf
  best_t1 <- NA_real_
  best_t2 <- NA_real_
  
  for (i1 in seq_along(thr_grid)) {
    t1 <- thr_grid[i1]
    for (i2 in seq_along(thr_grid)) {
      t2 <- thr_grid[i2]
      if (t2 <= t1) { pb_thr$tick(); next }
      
      pred_train <- predict_3class_threshold(sc_ok, t1, t2)
      acc <- mean(pred_train == yy_ok, na.rm = TRUE)
      
      if (is.finite(acc) && acc > best_acc) {
        best_acc <- acc
        best_t1 <- t1
        best_t2 <- t2
      }
      pb_thr$tick()
    }
  }
  
  cat(sprintf("\nBest thresholds in fold %d: t1=%.4f, t2=%.4f (train acc=%.4f)\n",
              f, best_t1, best_t2, best_acc))
  
  # 5) Score TEST patients (ALL regions) and classify
  cat("\nScoring TEST patients and predicting...\n")
  
  pb_tes <- progress_bar$new(
    format = "Test scoring [:bar] :current/:total (:percent) | ETA: :eta",
    total  = length(test_pat),
    clear  = FALSE,
    width  = 70
  )
  
  score_test <- rep(NA_real_, length(test_pat))
  
  for (j in seq_along(test_pat)) {
    cnt <- count_patient_tiers(
      df = test_df, rid = test_pat[j],
      regions = region_names, grid = grid,
      env_list = env_list,
      a = a, b = b,
      max_visits = max_visits, use_pct_change = use_pct_change
    )
    score_test[j] <- score_from_counts(cnt["n_ring"], cnt["n_out"], cnt["n_tot"], best_w1, best_w2)
    pb_tes$tick()
  }
  
  pred_test <- predict_3class_threshold(score_test, best_t1, best_t2)
  
  # fallback for missing scores: majority in TRAIN fold
  maj_train <- as.numeric(names(sort(table(y_train), decreasing = TRUE))[1])
  pred_test[!is.finite(score_test)] <- maj_train
  
  all_true <- c(all_true, y_test)
  all_pred <- c(all_pred, pred_test)
  
  cat(sprintf("\nFold %d done (end: %s)\n", f, format(Sys.time(), "%H:%M:%S")))
}

cat(sprintf("\nTotal CV time: %.1f minutes\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

cat("\n==================== OVERALL RESULTS (CV aggregated) ====================\n")

# ============================================================
# Functional CN-envelope ordinal deviation classifier (Term 1)
# 3-class: CN vs MCI vs AD
# Patient-level stratified 5-fold CV
#
# Core idea:
#  - Build CN median/MAD envelopes over time (per region, per time grid point)
#  - For each patient, count how many (region,time) points fall:
#       CN band -> 0, Ring -> w1, Outside -> w2
#  - Score = (w1*n_ring + w2*n_out) / n_total
#  - Tune (w1,w2) inside each fold (FAST: precompute counts once)
#  - Tune 2 thresholds (t1,t2) on TRAIN for max accuracy
#
# Progress bars:
#  - Envelope building (per region)
#  - Tier-count precompute (per patient)
#  - (w1,w2) search (coarse + refine)
#  - Threshold search (t1,t2)
#  - Test scoring (per patient)
#
# Notes:
#  - Uses <=4 earliest visits per patient
#  - Uses percent-change normalization per patient by default
#  - Optionally subsamples regions for tuning to reduce load
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(progress)
})

set.seed(123)

# ----------------------------
# USER SETTINGS
# ----------------------------
K_FOLDS        <- 5
max_visits     <- 4
region_start   <- 17
classes        <- c(1, 2, 4)      # CN, MCI, AD
grid_step      <- 6              # months
Tmax_use       <- 48             # cap horizon
use_pct_change <- TRUE           # recommended

# Envelope widths (fixed)
a <- 1.0                          # inner CN band multiplier
b <- 2.5                          # outer band multiplier (ring is between a and b)

# -------- Penalty tuning grids (FAST) --------
# Coarse search (small)
w1_vals_coarse <- seq(0, 10, by = 1.0)
w2_vals_coarse <- seq(0, 15, by = 1.0)

# Refine around best (still fast)
w1_refine_window <- 1.5
w2_refine_window <- 2.0
w1_refine_step   <- 0.25
w2_refine_step   <- 0.25

# Threshold tuning grid size (quantile grid)
n_thr_grid <- 80

# -------- Optional: reduce regions for tuning only (keeps final scoring on ALL regions) --------
# If you have tons of regions, tuning can still be slower due to tier-count precompute.
# This reduces that cost while keeping final chosen (w1,w2) applied to all regions.
use_region_subset_for_tuning <- TRUE
regions_tune_max <- 40   # e.g. 30-60 is reasonable

# ----------------------------
# 0) Load + basic cleanup
# ----------------------------
dat <- read.csv("data.csv")
dat$DIAGNOSIS <- as.numeric(as.character(dat$DIAGNOSIS))
dat$month     <- as.numeric(as.character(dat$month))

region_cols  <- region_start:ncol(dat)
region_names <- colnames(dat)[region_cols]

df0 <- dat %>%
  filter(!is.na(RID), !is.na(month), !is.na(DIAGNOSIS)) %>%
  filter(DIAGNOSIS %in% classes) %>%
  select(RID, month, DIAGNOSIS, all_of(region_names))

# Ensure unique names
colnames(df0) <- make.unique(colnames(df0))
region_names  <- colnames(df0)[4:ncol(df0)]

label_class <- function(x) factor(x, levels = c(1,2,4), labels = c("CN","MCI","AD"))

# ----------------------------
# Helper: patient-level stratified folds
# ----------------------------
patient_labels <- df0 %>% distinct(RID, DIAGNOSIS)

make_folds <- function(df_rid_lab, k = 5) {
  folds <- vector("list", k)
  for (i in 1:k) folds[[i]] <- integer(0)
  
  for (cl in sort(unique(df_rid_lab$DIAGNOSIS))) {
    rids <- df_rid_lab$RID[df_rid_lab$DIAGNOSIS == cl]
    rids <- sample(rids)
    parts <- split(rids, rep(1:k, length.out = length(rids)))
    for (i in 1:k) folds[[i]] <- c(folds[[i]], parts[[i]])
  }
  folds
}

folds <- make_folds(patient_labels, K_FOLDS)

# ----------------------------
# Helper: interpolation to grid
# ----------------------------
interp_to_grid <- function(t, y, grid) {
  if (length(t) < 2) return(rep(NA_real_, length(grid)))
  approx(x = t, y = y, xout = grid, method = "linear", rule = 1, ties = "ordered")$y
}

# ----------------------------
# Helper: build patient curves matrix for one region (n_pat x n_grid)
# using <= max_visits visits per patient
# ----------------------------
build_region_curves <- function(df, patient_ids, region, grid,
                                max_visits = 4, use_pct_change = TRUE) {
  
  df_sub <- df %>%
    filter(RID %in% patient_ids) %>%
    select(RID, month, all_of(region)) %>%
    rename(y = all_of(region)) %>%
    filter(!is.na(y)) %>%
    arrange(RID, month) %>%
    group_by(RID) %>%
    slice_head(n = max_visits) %>%
    ungroup()
  
  out <- matrix(NA_real_, nrow = length(patient_ids), ncol = length(grid))
  rownames(out) <- as.character(patient_ids)
  
  split_pat <- split(df_sub, df_sub$RID)
  
  for (i in seq_along(patient_ids)) {
    rid <- patient_ids[i]
    si <- split_pat[[as.character(rid)]]
    if (is.null(si) || nrow(si) < 2) next
    
    t <- si$month
    y <- si$y
    
    if (use_pct_change) {
      if (!is.finite(y[1]) || y[1] == 0) next
      y <- 100 * (y - y[1]) / y[1]
    }
    
    out[i, ] <- interp_to_grid(t, y, grid)
  }
  
  out
}

# ----------------------------
# Helper: build CN envelopes for one region
# Returns list(m = median_t, s = MAD_t) each length n_grid
# ----------------------------
build_cn_envelope_region <- function(curves_cn_mat) {
  m <- apply(curves_cn_mat, 2, median, na.rm = TRUE)
  s <- apply(curves_cn_mat, 2, mad,    na.rm = TRUE)
  s[!is.finite(s) | s == 0] <- 1e-6
  list(m = m, s = s)
}

# ----------------------------
# Helper: for ONE patient, compute counts of points in CN / ring / outside
# across ALL regions and time grid (using env_list)
# Returns c(n_cn, n_ring, n_out, n_tot)
# ----------------------------
count_patient_tiers <- function(df, rid, regions, grid, env_list,
                                a, b,
                                max_visits = 4, use_pct_change = TRUE) {
  
  subp <- df %>%
    filter(RID == rid) %>%
    arrange(month) %>%
    slice_head(n = max_visits)
  
  if (nrow(subp) < 2) return(c(n_cn=0L, n_ring=0L, n_out=0L, n_tot=0L))
  
  n_cn <- 0L; n_ring <- 0L; n_out <- 0L; n_tot <- 0L
  tvec_all <- subp$month
  
  for (reg in regions) {
    yvec <- subp[[reg]]
    ok <- which(!is.na(yvec))
    if (length(ok) < 2) next
    
    t <- tvec_all[ok]
    y <- yvec[ok]
    
    if (use_pct_change) {
      if (!is.finite(y[1]) || y[1] == 0) next
      y <- 100 * (y - y[1]) / y[1]
    }
    
    y_grid <- interp_to_grid(t, y, grid)
    good <- is.finite(y_grid)
    if (!any(good)) next
    
    env <- env_list[[reg]]
    m <- env$m; s <- env$s
    
    cn_low   <- m - a * s
    cn_high  <- m + a * s
    out_low  <- m - b * s
    out_high <- m + b * s
    
    inside_cn  <- good & (y_grid >= cn_low)  & (y_grid <= cn_high)
    inside_out <- good & (y_grid >= out_low) & (y_grid <= out_high)
    
    ring <- inside_out & !inside_cn
    out  <- good & !inside_out
    
    n_cn   <- n_cn   + sum(inside_cn)
    n_ring <- n_ring + sum(ring)
    n_out  <- n_out  + sum(out)
    n_tot  <- n_tot  + sum(good)
  }
  
  c(n_cn = n_cn, n_ring = n_ring, n_out = n_out, n_tot = n_tot)
}

# ----------------------------
# Helper: score from counts (vectorized-friendly)
# ----------------------------
score_from_counts <- function(n_ring, n_out, n_tot, w1, w2) {
  sc <- (w1 * n_ring + w2 * n_out) / n_tot
  sc[!is.finite(sc)] <- NA_real_
  sc
}

# ----------------------------
# Helper: thresholding for 3-class from a 1D score
# pred = CN if score <= t1; MCI if t1 < score <= t2; AD if score > t2
# ----------------------------
predict_3class_threshold <- function(score, t1, t2) {
  ifelse(score <= t1, 1,
         ifelse(score <= t2, 2, 4))
}

# ============================================================
# CV loop
# ============================================================
all_true <- integer(0)
all_pred <- integer(0)

cat("Running functional CN-envelope Term1 classifier (3-class) ...\n")
t0 <- Sys.time()

for (f in 1:K_FOLDS) {
  
  cat("\n============================================================\n")
  cat(sprintf("FOLD %d/%d (start: %s)\n", f, K_FOLDS, format(Sys.time(), "%H:%M:%S")))
  cat("============================================================\n")
  
  test_rids  <- folds[[f]]
  train_rids <- setdiff(patient_labels$RID, test_rids)
  
  train_df <- df0 %>% filter(RID %in% train_rids)
  test_df  <- df0 %>% filter(RID %in% test_rids)
  
  train_pat <- sort(unique(train_df$RID))
  test_pat  <- sort(unique(test_df$RID))
  
  # time grid from training only
  train_first4 <- train_df %>%
    arrange(RID, month) %>%
    group_by(RID) %>%
    slice_head(n = max_visits) %>%
    ungroup()
  
  Tmax <- min(Tmax_use, max(train_first4$month, na.rm = TRUE))
  grid <- seq(0, Tmax, by = grid_step)
  cat("Grid: 0 to", Tmax, "by", grid_step, "months (n =", length(grid), ")\n")
  
  # CN patient ids in training
  cn_train_pat <- train_df %>% filter(DIAGNOSIS == 1) %>% distinct(RID) %>% pull(RID)
  if (length(cn_train_pat) < 5) {
    stop("Too few CN patients in training fold to build functional envelopes robustly.")
  }
  
  # Labels per patient (TRAIN + TEST)
  y_train <- train_df %>%
    distinct(RID, DIAGNOSIS) %>%
    arrange(match(RID, train_pat)) %>%
    pull(DIAGNOSIS)
  
  y_test <- test_df %>%
    distinct(RID, DIAGNOSIS) %>%
    arrange(match(RID, test_pat)) %>%
    pull(DIAGNOSIS)
  
  # 1) Build CN envelopes per region (TRAIN only)
  cat("Building CN envelopes per region (TRAIN only)...\n")
  
  pb_env <- progress_bar$new(
    format = "Envelopes [:bar] :current/:total (:percent) | ETA: :eta",
    total  = length(region_names),
    clear  = FALSE,
    width  = 70
  )
  
  env_list <- vector("list", length(region_names))
  names(env_list) <- region_names
  
  for (reg in region_names) {
    cn_curves <- build_region_curves(train_df, cn_train_pat, reg, grid,
                                     max_visits = max_visits,
                                     use_pct_change = use_pct_change)
    env_list[[reg]] <- build_cn_envelope_region(cn_curves)
    pb_env$tick()
  }
  
  # Choose regions used for tuning (optional)
  if (use_region_subset_for_tuning && length(region_names) > regions_tune_max) {
    set.seed(1000 + f)  # fold-stable
    regions_tune <- sample(region_names, regions_tune_max)
    cat(sprintf("\nTuning penalties using a subset of %d/%d regions.\n",
                length(regions_tune), length(region_names)))
  } else {
    regions_tune <- region_names
    cat(sprintf("\nTuning penalties using ALL %d regions.\n", length(region_names)))
  }
  
  # 2) FAST tuning (w1,w2): precompute tier counts once, then cheap grid search
  cat("\nPrecomputing tier counts for TRAIN patients (once)...\n")
  pb_cnt <- progress_bar$new(
    format = "Tier counts [:bar] :current/:total (:percent) | ETA: :eta",
    total  = length(train_pat),
    clear  = FALSE,
    width  = 70
  )
  
  cnt_mat <- matrix(0, nrow = length(train_pat), ncol = 4)
  colnames(cnt_mat) <- c("n_cn","n_ring","n_out","n_tot")
  
  for (j in seq_along(train_pat)) {
    cnt_mat[j, ] <- count_patient_tiers(
      df = train_df, rid = train_pat[j],
      regions = regions_tune, grid = grid,
      env_list = env_list,
      a = a, b = b,
      max_visits = max_visits, use_pct_change = use_pct_change
    )
    pb_cnt$tick()
  }
  
  n_ring <- cnt_mat[, "n_ring"]
  n_out  <- cnt_mat[, "n_out"]
  n_tot  <- cnt_mat[, "n_tot"]
  
  ok_pat <- is.finite(y_train) & (n_tot > 0)
  if (sum(ok_pat) < 10) stop("Too few TRAIN patients with valid tier counts for tuning.")
  
  yy <- y_train[ok_pat]
  nr <- n_ring[ok_pat]
  no <- n_out[ok_pat]
  nt <- n_tot[ok_pat]
  
  # ---- Coarse grid search (w2 > w1 only) ----
  cat("\nTuning penalties (coarse grid, w2>w1)...\n")
  grid_w <- expand.grid(w1 = w1_vals_coarse, w2 = w2_vals_coarse) %>%
    filter(w2 > w1)
  
  pb_w <- progress_bar$new(
    format = "Penalty search (coarse) [:bar] :current/:total (:percent) | ETA: :eta",
    total  = nrow(grid_w),
    clear  = FALSE,
    width  = 70
  )
  
  best_sep <- -Inf
  best_w1  <- NA_real_
  best_w2  <- NA_real_
  
  for (i in seq_len(nrow(grid_w))) {
    w1 <- grid_w$w1[i]
    w2 <- grid_w$w2[i]
    
    sc <- (w1 * nr + w2 * no) / nt
    
    group_means  <- tapply(sc, yy, mean, na.rm = TRUE)
    overall_mean <- mean(sc, na.rm = TRUE)
    n_g <- table(yy)
    
    between_var <- sum(n_g * (group_means - overall_mean)^2, na.rm = TRUE)
    within_var  <- sum(tapply(sc, yy, var, na.rm = TRUE), na.rm = TRUE)
    
    sep <- between_var / (within_var + 1e-12)
    
    if (is.finite(sep) && sep > best_sep) {
      best_sep <- sep
      best_w1  <- w1
      best_w2  <- w2
    }
    pb_w$tick()
  }
  
  cat(sprintf("\nBest coarse penalties: w1=%.3f, w2=%.3f (sep=%.4f)\n",
              best_w1, best_w2, best_sep))
  
  # ---- Refine search around best ----
  cat("\nRefining penalties around best (still w2>w1)...\n")
  w1_lo <- max(0,  best_w1 - w1_refine_window)
  w1_hi <- min(max(w1_vals_coarse), best_w1 + w1_refine_window)
  w2_lo <- max(0,  best_w2 - w2_refine_window)
  w2_hi <- min(max(w2_vals_coarse), best_w2 + w2_refine_window)
  
  w1_vals_fine <- seq(w1_lo, w1_hi, by = w1_refine_step)
  w2_vals_fine <- seq(w2_lo, w2_hi, by = w2_refine_step)
  
  grid_w2 <- expand.grid(w1 = w1_vals_fine, w2 = w2_vals_fine) %>%
    filter(w2 > w1)
  
  pb_w2 <- progress_bar$new(
    format = "Penalty search (refine) [:bar] :current/:total (:percent) | ETA: :eta",
    total  = nrow(grid_w2),
    clear  = FALSE,
    width  = 70
  )
  
  best_sep2 <- -Inf
  best_w1_2 <- best_w1
  best_w2_2 <- best_w2
  
  for (i in seq_len(nrow(grid_w2))) {
    w1 <- grid_w2$w1[i]
    w2 <- grid_w2$w2[i]
    
    sc <- (w1 * nr + w2 * no) / nt
    
    group_means  <- tapply(sc, yy, mean, na.rm = TRUE)
    overall_mean <- mean(sc, na.rm = TRUE)
    n_g <- table(yy)
    
    between_var <- sum(n_g * (group_means - overall_mean)^2, na.rm = TRUE)
    within_var  <- sum(tapply(sc, yy, var, na.rm = TRUE), na.rm = TRUE)
    
    sep <- between_var / (within_var + 1e-12)
    
    if (is.finite(sep) && sep > best_sep2) {
      best_sep2 <- sep
      best_w1_2 <- w1
      best_w2_2 <- w2
    }
    pb_w2$tick()
  }
  
  best_w1 <- best_w1_2
  best_w2 <- best_w2_2
  best_sep <- best_sep2
  
  cat(sprintf("\nBest penalties in fold %d: w1=%.3f, w2=%.3f (sep=%.4f)\n",
              f, best_w1, best_w2, best_sep))
  
  # 3) Score TRAIN patients with best (w1,w2) using ALL regions (fresh counts)
  cat("\nScoring TRAIN patients with best (w1,w2) using ALL regions...\n")
  pb_trs <- progress_bar$new(
    format = "Train scoring [:bar] :current/:total (:percent) | ETA: :eta",
    total  = length(train_pat),
    clear  = FALSE,
    width  = 70
  )
  
  score_train <- rep(NA_real_, length(train_pat))
  
  for (j in seq_along(train_pat)) {
    cnt <- count_patient_tiers(
      df = train_df, rid = train_pat[j],
      regions = region_names, grid = grid,
      env_list = env_list,
      a = a, b = b,
      max_visits = max_visits, use_pct_change = use_pct_change
    )
    score_train[j] <- score_from_counts(cnt["n_ring"], cnt["n_out"], cnt["n_tot"], best_w1, best_w2)
    pb_trs$tick()
  }
  
  # 4) Tune thresholds (t1,t2) on TRAIN to maximize accuracy
  cat("\nTuning thresholds (t1,t2) on TRAIN for max accuracy...\n")
  
  ok <- is.finite(score_train)
  sc_ok <- score_train[ok]
  yy_ok <- y_train[ok]
  
  qs <- seq(0.02, 0.98, length.out = n_thr_grid)
  thr_grid <- as.numeric(quantile(sc_ok, probs = qs, na.rm = TRUE, names = FALSE))
  thr_grid <- sort(unique(thr_grid))
  
  if (length(thr_grid) < 5) stop("Too few distinct training scores to tune thresholds.")
  
  pb_thr <- progress_bar$new(
    format = "Threshold search [:bar] :current/:total (:percent) | ETA: :eta",
    total  = length(thr_grid) * length(thr_grid),
    clear  = FALSE,
    width  = 70
  )
  
  best_acc <- -Inf
  best_t1 <- NA_real_
  best_t2 <- NA_real_
  
  for (i1 in seq_along(thr_grid)) {
    t1 <- thr_grid[i1]
    for (i2 in seq_along(thr_grid)) {
      t2 <- thr_grid[i2]
      if (t2 <= t1) { pb_thr$tick(); next }
      
      pred_train <- predict_3class_threshold(sc_ok, t1, t2)
      acc <- mean(pred_train == yy_ok, na.rm = TRUE)
      
      if (is.finite(acc) && acc > best_acc) {
        best_acc <- acc
        best_t1 <- t1
        best_t2 <- t2
      }
      pb_thr$tick()
    }
  }
  
  cat(sprintf("\nBest thresholds in fold %d: t1=%.4f, t2=%.4f (train acc=%.4f)\n",
              f, best_t1, best_t2, best_acc))
  
  # 5) Score TEST patients (ALL regions) and classify
  cat("\nScoring TEST patients and predicting...\n")
  
  pb_tes <- progress_bar$new(
    format = "Test scoring [:bar] :current/:total (:percent) | ETA: :eta",
    total  = length(test_pat),
    clear  = FALSE,
    width  = 70
  )
  
  score_test <- rep(NA_real_, length(test_pat))
  
  for (j in seq_along(test_pat)) {
    cnt <- count_patient_tiers(
      df = test_df, rid = test_pat[j],
      regions = region_names, grid = grid,
      env_list = env_list,
      a = a, b = b,
      max_visits = max_visits, use_pct_change = use_pct_change
    )
    score_test[j] <- score_from_counts(cnt["n_ring"], cnt["n_out"], cnt["n_tot"], best_w1, best_w2)
    pb_tes$tick()
  }
  
  pred_test <- predict_3class_threshold(score_test, best_t1, best_t2)
  
  # fallback for missing scores: majority in TRAIN fold
  maj_train <- as.numeric(names(sort(table(y_train), decreasing = TRUE))[1])
  pred_test[!is.finite(score_test)] <- maj_train
  
  all_true <- c(all_true, y_test)
  all_pred <- c(all_pred, pred_test)
  
  cat(sprintf("\nFold %d done (end: %s)\n", f, format(Sys.time(), "%H:%M:%S")))
}

cat(sprintf("\nTotal CV time: %.1f minutes\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

cat("\n==================== OVERALL RESULTS (CV aggregated) ====================\n")
cm <- table(True = label_class(all_true), Pred = label_class(all_pred))
print(cm)

acc <- mean(all_true == all_pred)
cat(sprintf("\nOverall CV accuracy: %.4f\n", acc))
cat("=======================================================================\n")




