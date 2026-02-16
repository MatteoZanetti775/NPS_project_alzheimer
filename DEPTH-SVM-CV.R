#---- Libraries ----
library(ddalpha)
library(dplyr)
library(tidyr)
library(e1071)   # SVM
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

#---- 2) Robust standardization (median + MAD) ----
robust_center <- apply(data[feature_cols], 2, median, na.rm = TRUE)
robust_scale  <- apply(data[feature_cols], 2, mad,    na.rm = TRUE)
robust_scale[robust_scale == 0] <- 1e-6

data_std <- data
data_std[feature_cols] <- sweep(data_std[feature_cols], 2, robust_center, "-")
data_std[feature_cols] <- sweep(data_std[feature_cols], 2, robust_scale,  "/")

#---- 3) Compute projection depth space ----
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

depth_cols <- groups
print(table(depth_space$DIAGNOSIS))

#---- 4) Stratified K-fold helper (base R) ----
make_stratified_folds <- function(y, K = 5, seed = 123) {
  set.seed(seed)
  y <- droplevels(as.factor(y))
  folds <- vector("list", K)
  for (k in seq_len(K)) folds[[k]] <- integer(0)
  
  for (lev in levels(y)) {
    idx <- which(y == lev)
    idx <- sample(idx)
    parts <- split(idx, rep(1:K, length.out = length(idx)))
    for (k in seq_len(K)) folds[[k]] <- c(folds[[k]], parts[[k]])
  }
  
  lapply(folds, sort)
}

#---- 5) Cross-validated SVM in depth space ----
cv_svm_depth <- function(depth_df,
                         depth_cols,
                         label_col = "DIAGNOSIS",
                         K = 5,
                         cost = 1,
                         gamma = 0.5,
                         seed = 123) {
  
  folds <- make_stratified_folds(depth_df[[label_col]], K = K, seed = seed)
  
  acc_vec <- numeric(K)
  all_pred <- rep(NA, nrow(depth_df))
  all_true <- depth_df[[label_col]]
  
  for (fold_id in seq_len(K)) {
    test_idx  <- folds[[fold_id]]
    train_idx <- setdiff(seq_len(nrow(depth_df)), test_idx)
    
    trainD <- depth_df[train_idx, ]
    testD  <- depth_df[test_idx,  ]
    
    svm_model <- svm(
      DIAGNOSIS ~ .,
      data = trainD[, c(depth_cols, "DIAGNOSIS")],
      kernel = "radial",
      cost = cost,
      gamma = gamma
    )
    
    pred <- predict(svm_model, newdata = testD[, depth_cols])
    acc_vec[fold_id] <- mean(pred == testD[[label_col]])
    all_pred[test_idx] <- as.character(pred)
  }
  
  all_pred <- factor(all_pred, levels = levels(all_true))
  
  list(
    mean_accuracy = mean(acc_vec),
    sd_accuracy   = sd(acc_vec),
    confusion     = table(Predicted = all_pred, True = all_true)
  )
}

#---- 6) Run 5-fold CV SVM ----
K_folds <- 5

svm_cv_res <- cv_svm_depth(
  depth_df = depth_space,
  depth_cols = depth_cols,
  label_col = "DIAGNOSIS",
  K = K_folds,
  cost = 1,
  gamma = 0.5,
  seed = 123
)

cat("\nSVM 5-fold CV mean accuracy:", svm_cv_res$mean_accuracy, "\n")
cat("SVM 5-fold CV SD:", svm_cv_res$sd_accuracy, "\n\n")

cat("Cross-validated confusion matrix (SVM):\n")
print(svm_cv_res$confusion)

#---- 3D SVM decision boundary in depth space ----
library(rgl)
library(dplyr)

# Colors
cols_diag <- c("CN" = "blue", "MCI" = "orange", "AD" = "red")

# 1) Plot depth-space points
open3d()
points3d(
  depth_space$CN,
  depth_space$MCI,
  depth_space$AD,
  col   = cols_diag[depth_space$DIAGNOSIS],
  size  = 6,
  alpha = 0.6,
  xlab = "D_CN",
  ylab = "D_MCI",
  zlab = "D_AD"
)

# 2) Build a grid in depth space
res <- 40  # increase for smoother boundary if needed

cn_seq  <- seq(min(depth_space$CN),  max(depth_space$CN),  length.out = res)
mci_seq <- seq(min(depth_space$MCI), max(depth_space$MCI), length.out = res)
ad_seq  <- seq(min(depth_space$AD),  max(depth_space$AD),  length.out = res)

grid <- expand.grid(CN = cn_seq, MCI = mci_seq, AD = ad_seq)

#---- Fit SVM on full depth space (for visualization only) ----
svm_fit <- svm(
  DIAGNOSIS ~ .,
  data   = depth_space[, c("CN", "MCI", "AD", "DIAGNOSIS")],
  kernel = "radial",
  cost   = 1,
  gamma  = 0.5
)


# 3) Predict SVM class on grid
grid$pred <- predict(svm_fit, newdata = grid)

# 4) Extract approximate boundary points (class change)
boundary_pts <- grid %>%
  arrange(CN, MCI, AD) %>%
  group_by(CN, MCI) %>%
  mutate(bz = pred != lag(pred) | pred != lead(pred)) %>%
  ungroup() %>%
  arrange(MCI, AD, CN) %>%
  group_by(MCI, AD) %>%
  mutate(bx = pred != lag(pred) | pred != lead(pred)) %>%
  ungroup() %>%
  arrange(CN, AD, MCI) %>%
  group_by(CN, AD) %>%
  mutate(by = pred != lag(pred) | pred != lead(pred)) %>%
  ungroup() %>%
  filter(bx | by | bz)

# 5) Plot SVM decision boundary
points3d(
  boundary_pts$CN,
  boundary_pts$MCI,
  boundary_pts$AD,
  col   = "black",
  size  = 3,
  alpha = 0.4
)

#---- 7) Static 3D depth scatter + SVM decision boundary (export-friendly) ----
if (!requireNamespace("scatterplot3d", quietly = TRUE)) install.packages("scatterplot3d")
library(scatterplot3d)

# Colors wanted: AD red, MCI blue, CN green
pt_cols <- ifelse(depth_space$DIAGNOSIS == "AD",  "red",
                  ifelse(depth_space$DIAGNOSIS == "MCI", "blue",
                         ifelse(depth_space$DIAGNOSIS == "CN",  "green", "gray")))

# Fit SVM on full depth space (visualization only)
svm_fit <- svm(
  DIAGNOSIS ~ .,
  data   = depth_space[, c("CN", "MCI", "AD", "DIAGNOSIS")],
  kernel = "radial",
  cost   = 1,
  gamma  = 0.5
)

# Build a 3D grid in depth space
res <- 125  # increase for smoother boundary (slower + heavier)
cn_seq  <- seq(min(depth_space$CN),  max(depth_space$CN),  length.out = res)
mci_seq <- seq(min(depth_space$MCI), max(depth_space$MCI), length.out = res)
ad_seq  <- seq(min(depth_space$AD),  max(depth_space$AD),  length.out = res)

grid <- expand.grid(CN = cn_seq, MCI = mci_seq, AD = ad_seq)
grid$pred <- predict(svm_fit, newdata = grid)

# Convert predictions to an array: [CN, MCI, AD]
pred_arr <- array(as.integer(grid$pred), dim = c(res, res, res))

# Boundary mask: a grid point is boundary if it differs from any +1 neighbor
bmask <- array(FALSE, dim = c(res, res, res))
bmask[-res, , ] <- bmask[-res, , ] | (pred_arr[-res, , ] != pred_arr[-1, , ])
bmask[, -res, ] <- bmask[, -res, ] | (pred_arr[, -res, ] != pred_arr[, -1, ])
bmask[, , -res] <- bmask[, , -res] | (pred_arr[, , -res] != pred_arr[, , -1])

# Extract boundary coordinates
idx <- which(bmask, arr.ind = TRUE)
boundary_pts <- data.frame(
  CN  = cn_seq[idx[, 1]],
  MCI = mci_seq[idx[, 2]],
  AD  = ad_seq[idx[, 3]]
)

# --- Plot: points + boundary overlay (static) ---
s3d <- scatterplot3d(
  x = depth_space$AD, y = depth_space$MCI, z = depth_space$CN,
  color = pt_cols,
  pch = 16,
  cex.symbols = 0.4,   # make smaller/larger as you like
  xlab = "Depth wrt AD",
  ylab = "Depth wrt MCI",
  zlab = "Depth wrt CN",
  main = "Depth space + SVM decision boundary (approx.)"
)

# Overlay boundary points (semi-transparent black)
s3d$points3d(
  x = boundary_pts$AD, y = boundary_pts$MCI, z = boundary_pts$CN,
  pch = 16,
  cex = 0.25,
  col = rgb(0, 0, 0, 0.25)
)



