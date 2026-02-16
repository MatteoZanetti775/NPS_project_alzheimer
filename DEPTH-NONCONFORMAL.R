#---- Libraries ----
library(ddalpha)
library(dplyr)
library(tidyr)
library(tibble)

#---- 0) Load cleaned data ----
data <- read.csv("dataQR_clean.csv")
data$DIAGNOSIS <- factor(data$DIAGNOSIS)
groups <- levels(data$DIAGNOSIS)
print(table(data$DIAGNOSIS))

#---- 1) Identify feature columns ----
feature_cols <- setdiff(
  names(data),
  c("RID", "DIAGNOSIS", "GENOTYPE", "DIAGNOSIS_NUM")
)
data[feature_cols] <- lapply(data[feature_cols], as.numeric)

#---- 2) Robust standardization (median + MAD on pooled data) ----
robust_center <- apply(data[feature_cols], 2, median, na.rm = TRUE)
robust_scale  <- apply(data[feature_cols], 2, mad,    na.rm = TRUE)
robust_scale[robust_scale == 0] <- 1e-6

data_std <- data
data_std[feature_cols] <- sweep(data_std[feature_cols], 2, robust_center, "-")
data_std[feature_cols] <- sweep(data_std[feature_cols], 2, robust_scale,  "/")

#---- 3) Compute projection depth space (D_CN, D_MCI, D_AD) ----
compute_depth_matrix_projection <- function(df, feature_cols,
                                            label_col = "DIAGNOSIS",
                                            seed = 123) {
  groups <- levels(df[[label_col]])
  X_all  <- as.matrix(df[, feature_cols, drop = FALSE])
  
  depth_mat <- matrix(NA_real_, nrow = nrow(df), ncol = length(groups))
  colnames(depth_mat) <- groups
  
  for (j in seq_along(groups)) {
    g  <- groups[j]
    Xg <- as.matrix(df[df[[label_col]] == g, feature_cols, drop = FALSE])
    set.seed(seed)
    depth_mat[, j] <- depth.projection(x = X_all, data = Xg)
  }
  
  depth_df <- as.data.frame(depth_mat)
  depth_df[[label_col]] <- df[[label_col]]
  depth_df
}

depth_space <- compute_depth_matrix_projection(
  df = data_std,
  feature_cols = feature_cols,
  label_col = "DIAGNOSIS",
  seed = 123
)

depth_space$DIAGNOSIS <- droplevels(factor(depth_space$DIAGNOSIS))
depth_cols <- levels(depth_space$DIAGNOSIS)
depth_space[depth_cols] <- lapply(depth_space[depth_cols], as.numeric)
depth_space <- depth_space[complete.cases(depth_space[, c(depth_cols, "DIAGNOSIS")]), ]

cat("\nDepth space built. Class counts:\n")
print(table(depth_space$DIAGNOSIS))

#---- 4) Stratified train/calibration/test split (base R) ----
# Inductive conformal classification needs:
#   - proper training set (to compute class depths)
#   - calibration set (to compute quantiles of nonconformity)
#   - test set (to output conformal prediction sets)

stratified_split_threeway <- function(y, p_train = 0.6, p_cal = 0.2, seed = 123) {
  set.seed(seed)
  y <- droplevels(as.factor(y))
  
  idx_train <- integer(0)
  idx_cal   <- integer(0)
  idx_test  <- integer(0)
  
  for (lev in levels(y)) {
    idx <- which(y == lev)
    idx <- sample(idx)
    
    n <- length(idx)
    n_train <- max(1, floor(p_train * n))
    n_cal   <- max(1, floor(p_cal   * n))
    # remainder goes to test
    n_test  <- n - n_train - n_cal
    if (n_test < 1) {
      # ensure at least 1 test if possible by stealing from train
      if (n_train > 1) {
        n_train <- n_train - 1
        n_test  <- n_test + 1
      } else if (n_cal > 1) {
        n_cal <- n_cal - 1
        n_test <- n_test + 1
      }
    }
    
    idx_train <- c(idx_train, idx[1:n_train])
    idx_cal   <- c(idx_cal,   idx[(n_train + 1):(n_train + n_cal)])
    idx_test  <- c(idx_test,  idx[(n_train + n_cal + 1):n])
  }
  
  list(
    train = sort(idx_train),
    cal   = sort(idx_cal),
    test  = sort(idx_test)
  )
}

splits <- stratified_split_threeway(depth_space$DIAGNOSIS, p_train = 0.6, p_cal = 0.2, seed = 123)

trainD <- depth_space[splits$train, , drop = FALSE]
calD   <- depth_space[splits$cal,   , drop = FALSE]
testD  <- depth_space[splits$test,  , drop = FALSE]

trainD$DIAGNOSIS <- droplevels(trainD$DIAGNOSIS)
calD$DIAGNOSIS   <- factor(calD$DIAGNOSIS,  levels = levels(trainD$DIAGNOSIS))
testD$DIAGNOSIS  <- factor(testD$DIAGNOSIS, levels = levels(trainD$DIAGNOSIS))

cat("\nSplit sizes:\n")
cat("  Train:", nrow(trainD), "\n")
cat("  Cal  :", nrow(calD), "\n")
cat("  Test :", nrow(testD), "\n")

cat("\nClass counts (train/cal/test):\n")
print(table(trainD$DIAGNOSIS))
print(table(calD$DIAGNOSIS))
print(table(testD$DIAGNOSIS))

#---- 5) Option A nonconformity: A(x,y) = - D_y(x) ----
# Here D_y(x) is the projection depth w.r.t. the training data of class y.
# We'll compute depths of calibration+test points relative to each class in TRAIN.

compute_depths_wrt_train_classes <- function(train_std_df, new_std_df, feature_cols, label_col = "DIAGNOSIS", seed = 123) {
  groups <- levels(train_std_df[[label_col]])
  X_new  <- as.matrix(new_std_df[, feature_cols, drop = FALSE])
  
  depth_mat <- matrix(NA_real_, nrow = nrow(new_std_df), ncol = length(groups))
  colnames(depth_mat) <- groups
  
  for (j in seq_along(groups)) {
    g  <- groups[j]
    Xg <- as.matrix(train_std_df[train_std_df[[label_col]] == g, feature_cols, drop = FALSE])
    
    set.seed(seed)
    depth_mat[, j] <- depth.projection(x = X_new, data = Xg)
  }
  
  as.data.frame(depth_mat)
}

# IMPORTANT:
# The depth_space was computed using depths w.r.t. ALL data earlier.
# For conformal validity, we re-compute depths for cal/test using ONLY TRAIN as reference.

# We need the standardized ORIGINAL feature space split too:
data_std$DIAGNOSIS <- factor(data_std$DIAGNOSIS, levels = levels(depth_space$DIAGNOSIS))

train_std <- data_std[splits$train, , drop = FALSE]
cal_std   <- data_std[splits$cal,   , drop = FALSE]
test_std  <- data_std[splits$test,  , drop = FALSE]

train_std$DIAGNOSIS <- droplevels(train_std$DIAGNOSIS)
cal_std$DIAGNOSIS   <- factor(cal_std$DIAGNOSIS,  levels = levels(train_std$DIAGNOSIS))
test_std$DIAGNOSIS  <- factor(test_std$DIAGNOSIS, levels = levels(train_std$DIAGNOSIS))

# Depths of calibration points w.r.t. each TRAIN class
cal_depths <- compute_depths_wrt_train_classes(
  train_std_df = train_std,
  new_std_df   = cal_std,
  feature_cols = feature_cols,
  label_col    = "DIAGNOSIS",
  seed         = 123
)

# Depths of test points w.r.t. each TRAIN class
test_depths <- compute_depths_wrt_train_classes(
  train_std_df = train_std,
  new_std_df   = test_std,
  feature_cols = feature_cols,
  label_col    = "DIAGNOSIS",
  seed         = 123
)

#---- 6) Inductive conformal classification (per-class calibration) ----
# For each class y:
#   calibration nonconformity scores: A_i = - D_y(x_i) for calibration points with true label y
# For a new test point x:
#   p_y(x) = ( #{i in cal_y: A_i >= A_test(y)} + 1 ) / (n_cal_y + 1)
# Prediction set at level alpha: { y : p_y(x) > alpha }

conformal_p_values <- function(cal_depths, cal_labels, test_depths, alpha = 0.1) {
  classes <- levels(cal_labels)
  
  # Store calibration scores by class
  cal_scores <- vector("list", length(classes))
  names(cal_scores) <- classes
  
  for (y in classes) {
    idx_y <- which(cal_labels == y)
    # Option A: A = - depth_y
    cal_scores[[y]] <- -cal_depths[idx_y, y]
  }
  
  # Compute p-values for each test point and each class
  p_mat <- matrix(NA_real_, nrow = nrow(test_depths), ncol = length(classes))
  colnames(p_mat) <- classes
  
  for (y in classes) {
    scores_y <- cal_scores[[y]]
    n_y <- length(scores_y)
    if (n_y == 0) stop(paste("No calibration points for class", y))
    
    # test nonconformity under class y
    a_test <- -test_depths[, y]
    
    # p-value: (count >= a_test + 1) / (n_y + 1)
    # vectorized:
    p_mat[, y] <- sapply(a_test, function(at) (sum(scores_y >= at) + 1) / (n_y + 1))
  }
  
  # Prediction set: include labels with p > alpha
  pred_set <- apply(p_mat, 1, function(pv) names(pv)[pv > alpha])
  
  list(
    p_values = as.data.frame(p_mat),
    pred_set = pred_set
  )
}

# Choose miscoverage alpha (e.g., 0.1 -> 90% marginal coverage)
alpha <- 0.3

conf_res <- conformal_p_values(
  cal_depths  = cal_depths,
  cal_labels  = cal_std$DIAGNOSIS,
  test_depths = test_depths,
  alpha       = alpha
)

pvals_df <- conf_res$p_values
pred_sets <- conf_res$pred_set

#---- 7) Evaluate conformal classification: coverage and set sizes ----
true_test <- test_std$DIAGNOSIS

# Coverage: true label belongs to prediction set
covered <- mapply(function(S, y) y %in% S, pred_sets, as.character(true_test))
coverage <- mean(covered)

# Average size of prediction set
set_sizes <- sapply(pred_sets, length)
avg_size <- mean(set_sizes)

cat("\nConformal classification results (alpha =", alpha, ")\n")
cat("  Empirical coverage :", coverage, "\n")
cat("  Avg set size       :", avg_size, "\n")

cat("\nSet size distribution:\n")
print(table(set_sizes))

# Optional: how often we output a single label (most useful cases)
singleton_rate <- mean(set_sizes == 1)
cat("\nSingleton rate (|set|=1):", singleton_rate, "\n")

#---- 8) Also report 'point prediction' from conformal sets (optional) ----
# If singleton -> use that label; else NA (abstain)
point_pred <- sapply(pred_sets, function(S) if (length(S) == 1) S else NA_character_)
point_pred <- factor(point_pred, levels = levels(true_test))

usable <- !is.na(point_pred)
cat("\nAmong singleton predictions only:\n")
cat("  n singleton:", sum(usable), "out of", length(point_pred), "\n")
if (sum(usable) > 0) {
  conf_singleton <- table(Predicted = point_pred[usable], True = true_test[usable])
  print(conf_singleton)
  acc_singleton <- mean(point_pred[usable] == true_test[usable])
  cat("  Accuracy on singletons:", acc_singleton, "\n")
} else {
  cat("  (No singleton predictions at this alpha.)\n")
}

#---- 9) Inspect a few test examples (p-values + prediction sets) ----
# Combine into one table for quick look
inspect_df <- tibble::tibble(
  true_label = as.character(true_test),
  set_size   = set_sizes,
  pred_set   = sapply(pred_sets, function(S) paste(S, collapse = ", "))
) %>%
  dplyr::bind_cols(pvals_df)

print(head(inspect_df, 15))
