#---- Libraries ----
library(ddalpha)
library(dplyr)
library(tidyr)
library(class)    # knn
library(tibble)
library(ggplot2)

#---- 0) Load cleaned data ----
data <- read.csv("dataQR_clean.csv")

# Make sure DIAGNOSIS is a factor
data$DIAGNOSIS <- factor(data$DIAGNOSIS)
groups <- levels(data$DIAGNOSIS)
print(groups)
print(table(data$DIAGNOSIS))

#---- 1) Identify feature columns (all brain volumes) ----
feature_cols <- setdiff(
  names(data),
  c("RID", "DIAGNOSIS", "GENOTYPE", "DIAGNOSIS_NUM")
)

# Make sure all features are numeric
data[feature_cols] <- lapply(data[feature_cols], as.numeric)

# Optional quick sanity check
str(data[, c("DIAGNOSIS", feature_cols[1:3])])

#---- 2) Robust standardization (median + MAD on pooled data) ----
robust_center <- apply(data[feature_cols], 2, median, na.rm = TRUE)
robust_scale  <- apply(data[feature_cols], 2, mad,    na.rm = TRUE)
robust_scale[robust_scale == 0] <- 1e-6  # avoid division by zero

data_std <- data
data_std[feature_cols] <- sweep(data_std[feature_cols], 2, robust_center, "-")
data_std[feature_cols] <- sweep(data_std[feature_cols], 2, robust_scale,  "/")

#---- 3) Build depth space (projection depth per class) ----
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
    
    # ddalpha uses RNG internally for directions; fix it for reproducibility
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

print(head(depth_space))
print(table(depth_space$DIAGNOSIS))

depth_cols <- groups  # columns in depth_space corresponding to classes

#---- 4) Stratified K-fold CV helper (no extra packages) ----
make_stratified_folds <- function(y, K = 5, seed = 123) {
  set.seed(seed)
  y <- droplevels(as.factor(y))
  folds <- vector("list", K)
  for (k in seq_len(K)) folds[[k]] <- integer(0)
  
  for (lev in levels(y)) {
    idx <- which(y == lev)
    idx <- sample(idx, length(idx))          # shuffle within class
    parts <- split(idx, rep(1:K, length.out = length(idx)))
    for (k in seq_len(K)) folds[[k]] <- c(folds[[k]], parts[[k]])
  }
  
  lapply(folds, sort)
}

#---- 5) Cross-validated kNN in depth space (evaluate many k) ----
cv_knn_depth <- function(depth_df,
                         depth_cols,
                         label_col = "DIAGNOSIS",
                         k_values = seq(1, 101, by = 2),
                         K = 5,
                         seed = 123) {
  
  folds <- make_stratified_folds(depth_df[[label_col]], K = K, seed = seed)
  
  # Storage: rows = k values, cols = folds
  acc_mat <- matrix(NA_real_, nrow = length(k_values), ncol = K)
  colnames(acc_mat) <- paste0("fold", 1:K)
  rownames(acc_mat) <- paste0("k=", k_values)
  
  for (fold_id in seq_len(K)) {
    test_idx <- folds[[fold_id]]
    train_idx <- setdiff(seq_len(nrow(depth_df)), test_idx)
    
    trainD <- depth_df[train_idx, , drop = FALSE]
    testD  <- depth_df[test_idx,  , drop = FALSE]
    
    X_train <- trainD[, depth_cols, drop = FALSE]
    X_test  <- testD[,  depth_cols, drop = FALSE]
    y_train <- trainD[[label_col]]
    y_test  <- testD[[label_col]]
    
    for (i in seq_along(k_values)) {
      k <- k_values[i]
      
      pred <- knn(
        train = X_train,
        test  = X_test,
        cl    = y_train,
        k     = k
      )
      
      acc_mat[i, fold_id] <- mean(pred == y_test)
    }
  }
  
  acc_mean <- rowMeans(acc_mat, na.rm = TRUE)
  acc_sd   <- apply(acc_mat, 1, sd, na.rm = TRUE)
  
  data.frame(
    k = k_values,
    acc_mean = acc_mean,
    acc_sd = acc_sd
  )
}

# Choose K and k-grid
K_folds <- 5
k_values <- seq(1, 355, by = 2)   # odd k avoids ties

cv_results <- cv_knn_depth(
  depth_df = depth_space,
  depth_cols = depth_cols,
  label_col = "DIAGNOSIS",
  k_values = k_values,
  K = K_folds,
  seed = 123
)

print(head(cv_results))

# Best k by mean CV accuracy
best_row <- cv_results[which.max(cv_results$acc_mean), ]
cat("\nBest k (by mean CV accuracy):", best_row$k,
    "\nMean accuracy:", best_row$acc_mean,
    "\nSD across folds:", best_row$acc_sd, "\n")

#---- 6) Plot CV accuracy vs k ----
ggplot(cv_results, aes(x = k, y = acc_mean)) +
  geom_line(linewidth = 0.3) +
  geom_point(size = 1) +
  geom_ribbon(aes(ymin = acc_mean - acc_sd, ymax = acc_mean + acc_sd),
              alpha = 0.15, colour = NA) +
  labs(
    title = paste0("DD-kNN: ", K_folds, "-fold stratified CV accuracy vs k"),
    x = "Number of Neighbors (k)",
    y = "Mean CV Accuracy (± 1 SD)"
  ) +
  theme_bw() +
  theme(text = element_text(size = 14))

# ---- Overall average and SD of CV accuracy across k values ----
# ---- Accuracy summary for stable k region (k >= 80) ----
stable_k_min <- 1

cv_stable <- subset(cv_results, k >= stable_k_min)

avg_acc_stable <- mean(cv_stable$acc_mean)
sd_acc_stable  <- sd(cv_stable$acc_mean)

cat("\nCV accuracy for k >=", stable_k_min, ":\n")
cat("Average accuracy:", round(avg_acc_stable, 4), "\n")
cat("SD across k:",     round(sd_acc_stable, 4), "\n")
cat("Number of k values:", nrow(cv_stable), "\n")



#---- 7) Optional: confusion matrix at best k using CV folds (pooled preds) ----
# This gives you a single confusion matrix aggregated over all held-out folds.

cv_confusion_at_k <- function(depth_df, depth_cols, label_col = "DIAGNOSIS",
                              k = 5, K = 5, seed = 123) {
  
  folds <- make_stratified_folds(depth_df[[label_col]], K = K, seed = seed)
  
  all_pred <- rep(NA, nrow(depth_df))
  all_true <- depth_df[[label_col]]
  
  for (fold_id in seq_len(K)) {
    test_idx <- folds[[fold_id]]
    train_idx <- setdiff(seq_len(nrow(depth_df)), test_idx)
    
    trainD <- depth_df[train_idx, , drop = FALSE]
    testD  <- depth_df[test_idx,  , drop = FALSE]
    
    pred <- knn(
      train = trainD[, depth_cols, drop = FALSE],
      test  = testD[,  depth_cols, drop = FALSE],
      cl    = trainD[[label_col]],
      k     = k
    )
    
    all_pred[test_idx] <- as.character(pred)
  }
  
  all_pred <- factor(all_pred, levels = levels(depth_df[[label_col]]))
  table(Predicted = all_pred, True = all_true)
}

conf_mat_bestk <- cv_confusion_at_k(
  depth_df = depth_space,
  depth_cols = depth_cols,
  label_col = "DIAGNOSIS",
  k = best_row$k,
  K = K_folds,
  seed = 123
)

cat("\nCross-validated confusion matrix at best k:\n")
print(conf_mat_bestk)

#---- 8) Static 3D plot of depth space (AD red, MCI blue, CN green) ----
if (!requireNamespace("scatterplot3d", quietly = TRUE)) {
  install.packages("scatterplot3d")
}
library(scatterplot3d)

stopifnot(all(c("AD","MCI","CN") %in% colnames(depth_space)))

cols <- ifelse(depth_space$DIAGNOSIS == "AD",  "red",
               ifelse(depth_space$DIAGNOSIS == "MCI", "blue",
                      ifelse(depth_space$DIAGNOSIS == "CN",  "green", "gray")))

scatterplot3d(
  x = depth_space$AD,
  y = depth_space$MCI,
  z = depth_space$CN,
  color = cols,
  pch = 16,
  cex.symbols = 0.4,
  xlab = "Depth wrt AD",
  ylab = "Depth wrt MCI",
  zlab = "Depth wrt CN",
  main = "Depth space (Projection depth)"
)

legend(
  "topright",
  legend = c("AD", "MCI", "CN"),
  col = c("red", "blue", "green"),
  pch = 16,
  bty = "n"
)


#---- Static 3D depth scatter + kNN decision boundary (approx.) + progress bar ----
library(class)
if (!requireNamespace("scatterplot3d", quietly = TRUE)) install.packages("scatterplot3d")
if (!requireNamespace("progress", quietly = TRUE)) install.packages("progress")
library(scatterplot3d)
library(progress)

# Safety: ensure expected columns exist
stopifnot(all(c("CN","MCI","AD","DIAGNOSIS") %in% colnames(depth_space)))

# Choose k for kNN (use your CV best k if available, otherwise set manually)
k_knn <- 75

# Colors (as requested): AD red, MCI blue, CN green
pt_cols <- ifelse(depth_space$DIAGNOSIS == "AD",  "red",
                  ifelse(depth_space$DIAGNOSIS == "MCI", "blue",
                         ifelse(depth_space$DIAGNOSIS == "CN",  "green", "gray")))

# Training data in depth space
X_train <- depth_space[, c("CN","MCI","AD")]
y_train <- depth_space$DIAGNOSIS

# Build a 3D grid in depth space
res <- 75  # increase for smoother boundary (slower). Try 25 for speed.
cn_seq  <- seq(min(depth_space$CN),  max(depth_space$CN),  length.out = res)
mci_seq <- seq(min(depth_space$MCI), max(depth_space$MCI), length.out = res)
ad_seq  <- seq(min(depth_space$AD),  max(depth_space$AD),  length.out = res)

grid <- expand.grid(CN = cn_seq, MCI = mci_seq, AD = ad_seq)
X_test <- grid[, c("CN","MCI","AD")]
n_test <- nrow(X_test)

# ---- Predict kNN on grid WITH progress (chunked) ----
chunk_size <- 5000  # larger = faster, fewer updates; smaller = smoother bar
n_chunks <- ceiling(n_test / chunk_size)

pb <- progress_bar$new(
  format = "kNN grid prediction [:bar] :percent | ETA: :eta",
  total  = n_chunks,
  clear  = FALSE,
  width  = 60
)

pred_vec <- character(n_test)
starts <- seq(1, n_test, by = chunk_size)

for (s in starts) {
  e <- min(s + chunk_size - 1, n_test)
  
  pred_vec[s:e] <- as.character(knn(
    train = X_train,
    test  = X_test[s:e, , drop = FALSE],
    cl    = y_train,
    k     = k_knn
  ))
  
  pb$tick()
}

grid$pred <- factor(pred_vec, levels = levels(y_train))

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

# ---- Static plot: points + boundary overlay ----
s3d <- scatterplot3d(
  x = depth_space$AD, y = depth_space$MCI, z = depth_space$CN,
  color = pt_cols,
  pch = 16,
  cex.symbols = 0.4,  # smaller points; reduce further if needed
  xlab = "Depth wrt AD",
  ylab = "Depth wrt MCI",
  zlab = "Depth wrt CN",
  main = paste0("Depth space + kNN decision boundary (k=", k_knn, ", res=", res, ")")
)

# Overlay boundary points (semi-transparent black)
s3d$points3d(
  x = boundary_pts$AD, y = boundary_pts$MCI, z = boundary_pts$CN,
  pch = 16,
  cex = 0.25,
  col = rgb(0, 0, 0, 0.25)
)

legend(
  "topright",
  legend = c("AD", "MCI", "CN", "Boundary"),
  col = c("red", "blue", "green", rgb(0,0,0,0.6)),
  pch = 16,
  bty = "n"
)

