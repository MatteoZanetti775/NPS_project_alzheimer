library(ddalpha)   # already installed
library(dplyr)

data <- read.csv("dataQR_clean.csv")

## Make sure DIAGNOSIS is a factor with CN/MCI/AD
data$DIAGNOSIS <- factor(data$DIAGNOSIS)

groups <- levels(data$DIAGNOSIS)
groups


## ---------- Train / test split (e.g. 70% train, 30% test) ---------------

set.seed(123)  # for reproducibility

n <- nrow(data)
train_idx <- sample(seq_len(n), size = floor(0.7 * n))

train_data <- data[train_idx, , drop = FALSE]
test_data  <- data[-train_idx, , drop = FALSE]

table(train_data$DIAGNOSIS)
table(test_data$DIAGNOSIS)

## ---------- Spatial-depth classifier: max-depth rule ---------------------

# Use all brain volume columns, exclude ID / diagnosis / genotype
feature_cols <- setdiff(names(data), c("RID", "DIAGNOSIS", "GENOTYPE"))

depth_classify_spatial <- function(train_df, test_df, feature_cols, label_col = "DIAGNOSIS") {
  groups <- levels(train_df[[label_col]])
  X_test <- as.matrix(test_df[, feature_cols, drop = FALSE])
  
  depth_mat <- matrix(NA_real_,
                      nrow = nrow(test_df),
                      ncol = length(groups))
  colnames(depth_mat) <- groups
  
  for (j in seq_along(groups)) {
    g  <- groups[j]
    Xg <- as.matrix(train_df[train_df[[label_col]] == g, feature_cols, drop = FALSE])
    
    # depth of all test points w.r.t. group g
    depth_mat[, j] <- depth.spatial(x = X_test, data = Xg)
  }
  
  # assign each test point to group with maximum depth
  predicted <- groups[max.col(depth_mat, ties.method = "first")]
  
  list(
    depth_mat = depth_mat,
    predicted = factor(predicted, levels = groups)
  )
}

## Run classifier on our train/test split
clf_res <- depth_classify_spatial(train_data, test_data, feature_cols)

pred_labels <- clf_res$predicted
true_labels <- test_data$DIAGNOSIS

conf_mat <- table(Predicted = pred_labels, True = true_labels)
conf_mat

accuracy <- mean(pred_labels == true_labels)
accuracy

