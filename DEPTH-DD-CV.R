#---- Libraries ----
library(ddalpha)
library(dplyr)
library(tidyr)
library(rgl)

#---- 0) Load cleaned data ----
data <- read.csv("dataQR_top30.csv")
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

depth_space$DIAGNOSIS <- droplevels(factor(depth_space$DIAGNOSIS))
depth_cols <- levels(depth_space$DIAGNOSIS)   # should be c("CN","MCI","AD") (order depends on factor)
print(table(depth_space$DIAGNOSIS))

# Ensure numeric + complete cases
depth_space[depth_cols] <- lapply(depth_space[depth_cols], as.numeric)
depth_space <- depth_space[complete.cases(depth_space[, c(depth_cols, "DIAGNOSIS")]), ]

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

#---- 5) DD-alpha: K-fold CV on depth space ----
cv_ddalpha_depth <- function(depth_df,
                             depth_cols,
                             label_col = "DIAGNOSIS",
                             K = 5,
                             seed = 123) {
  
  folds <- make_stratified_folds(depth_df[[label_col]], K = K, seed = seed)
  
  acc_vec <- numeric(K)
  all_pred <- rep(NA, nrow(depth_df))
  all_true <- depth_df[[label_col]]
  
  for (fold_id in seq_len(K)) {
    test_idx  <- folds[[fold_id]]
    train_idx <- setdiff(seq_len(nrow(depth_df)), test_idx)
    
    trainD <- depth_df[train_idx, , drop = FALSE]
    testD  <- depth_df[test_idx,  , drop = FALSE]
    
    # ddalpha is picky: keep only present levels in training
    trainD[[label_col]] <- droplevels(trainD[[label_col]])
    testD[[label_col]]  <- factor(testD[[label_col]], levels = levels(trainD[[label_col]]))
    
    # Train on precomputed DD-plot
    train_ddplot <- trainD[, c(depth_cols, label_col), drop = FALSE]
    dd_model <- ddalpha.train(
      data      = train_ddplot,
      depth     = "ddplot",
      separator = "alpha"
    )
    
    # Predict held-out fold
    pred <- ddalpha.classify(dd_model, objects = testD[, depth_cols, drop = FALSE])
    pred <- factor(unlist(pred), levels = levels(trainD[[label_col]]))
    
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

#---- 6) Run 5-fold CV DD-alpha ----
K_folds <- 5

dd_cv_res <- cv_ddalpha_depth(
  depth_df = depth_space,
  depth_cols = depth_cols,
  label_col = "DIAGNOSIS",
  K = K_folds,
  seed = 123
)

cat("\nDD-alpha 5-fold CV mean accuracy:", dd_cv_res$mean_accuracy, "\n")
cat("DD-alpha 5-fold CV SD:", dd_cv_res$sd_accuracy, "\n\n")
cat("Cross-validated confusion matrix (DD-alpha):\n")
print(dd_cv_res$confusion)

#---- 7) 3D visualization: DD-alpha boundary (trained on full depth space) ----
# NOTE: This is for visualization only (CV results above are the performance estimate).

dd_model_full <- ddalpha.train(
  data      = depth_space[, c(depth_cols, "DIAGNOSIS"), drop = FALSE],
  depth     = "ddplot",
  separator = "alpha"
)

open3d()
cols_diag <- c("CN"="blue", "MCI"="orange", "AD"="red")

# Fallback if your factor levels aren't exactly CN/MCI/AD
cols_diag <- cols_diag[names(cols_diag) %in% levels(depth_space$DIAGNOSIS)]
if (length(cols_diag) == 0) {
  levs <- levels(depth_space$DIAGNOSIS)
  cols_diag <- setNames(rainbow(length(levs)), levs)
}

points3d(depth_space$CN, depth_space$MCI, depth_space$AD,
         col = cols_diag[depth_space$DIAGNOSIS], size = 6, alpha = 0.6,
         xlab = "D_CN", ylab = "D_MCI", zlab = "D_AD")

res <- 40
cn_seq  <- seq(min(depth_space$CN),  max(depth_space$CN),  length.out = res)
mci_seq <- seq(min(depth_space$MCI), max(depth_space$MCI), length.out = res)
ad_seq  <- seq(min(depth_space$AD),  max(depth_space$AD),  length.out = res)

grid <- expand.grid(CN = cn_seq, MCI = mci_seq, AD = ad_seq)

grid_pred <- ddalpha.classify(dd_model_full, objects = grid[, c("CN","MCI","AD")])
grid$pred <- factor(unlist(grid_pred), levels = levels(depth_space$DIAGNOSIS))

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

points3d(boundary_pts$CN, boundary_pts$MCI, boundary_pts$AD,
         col = "black", size = 2.5, alpha = 0.35)

centers <- depth_space %>%
  group_by(DIAGNOSIS) %>%
  summarize(CN = mean(CN), MCI = mean(MCI), AD = mean(AD), .groups = "drop")

for (i in seq_len(nrow(centers))) {
  text3d(centers$CN[i], centers$MCI[i], centers$AD[i],
         texts = as.character(centers$DIAGNOSIS[i]),
         color = cols_diag[as.character(centers$DIAGNOSIS[i])],
         cex = 2)
}
