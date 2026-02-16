#---- 0) Load cleaned data ----
data <- read.csv("dataQR_clean.csv")

## -----------------------------
## 0) Setup / sanity
## -----------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# Make sure DIAGNOSIS and GENOTYPE are treated as factors
data <- data %>%
  mutate(
    DIAGNOSIS = as.factor(DIAGNOSIS),
    GENOTYPE  = as.factor(GENOTYPE)
  )

# (Optional) If you want to force an order:
# data$DIAGNOSIS <- factor(data$DIAGNOSIS, levels = c("CN","MCI","AD"))

# Define feature columns: everything except identifiers/categoricals
id_cols <- c("RID", "DIAGNOSIS", "GENOTYPE")
feature_cols <- setdiff(names(data), id_cols)

# Ensure numeric features (just in case something was read as character)
# This will turn non-numeric to NA; you can inspect afterwards if needed.
data <- data %>%
  mutate(across(all_of(feature_cols), ~ suppressWarnings(as.numeric(.x))))


## -----------------------------
## 1) Kruskal-Wallis per region
## -----------------------------
kw_results <- map_dfr(feature_cols, function(feat) {
  df <- data %>% select(DIAGNOSIS, all_of(feat)) %>% drop_na()
  
  # Only run if all 3 classes are present after dropping NA
  if (n_distinct(df$DIAGNOSIS) < 3) {
    return(tibble(
      feature = feat,
      n = nrow(df),
      kw_statistic = NA_real_,
      kw_p = NA_real_
    ))
  }
  
  kt <- kruskal.test(df[[feat]] ~ df$DIAGNOSIS)
  tibble(
    feature = feat,
    n = nrow(df),
    kw_statistic = unname(kt$statistic),
    kw_p = kt$p.value
  )
}) %>%
  mutate(kw_p_adj_fdr = p.adjust(kw_p, method = "BH")) %>%
  arrange(kw_p_adj_fdr)

# Look at top results
print(head(kw_results, 20))


## -----------------------------
## 2) Post-hoc: pairwise Wilcoxon for significant regions
##    + flag regions significant in ALL 3 pairwise comparisons
## -----------------------------
alpha <- 0.05

sig_kw_features <- kw_results %>%
  filter(!is.na(kw_p_adj_fdr), kw_p_adj_fdr < alpha) %>%
  pull(feature)

pairwise_results <- map_dfr(sig_kw_features, function(feat) {
  df <- data %>% select(DIAGNOSIS, all_of(feat)) %>% drop_na()
  
  pw <- pairwise.wilcox.test(
    x = df[[feat]],
    g = df$DIAGNOSIS,
    p.adjust.method = "BH",
    exact = FALSE
  )
  
  # Convert the matrix of p-values to long form
  pmat <- pw$p.value
  if (is.null(pmat)) return(tibble())  # safety
  
  as.data.frame(as.table(pmat)) %>%
    rename(group1 = Var1, group2 = Var2, p_adj = Freq) %>%
    mutate(feature = feat) %>%
    filter(!is.na(p_adj))
})

# Determine if a feature is significant in all three pairwise comparisons
# Note: pairwise.wilcox.test returns only unique pairs; we count them.
all_pairs_needed <- 3  # CN-MCI, CN-AD, MCI-AD

feature_pairwise_summary <- pairwise_results %>%
  mutate(sig = p_adj < alpha) %>%
  group_by(feature) %>%
  summarise(
    n_pairs_tested = n(),
    n_pairs_sig = sum(sig),
    sig_all_pairs = (n_pairs_sig == all_pairs_needed),
    .groups = "drop"
  ) %>%
  arrange(desc(sig_all_pairs), desc(n_pairs_sig))

# Join summary back with KW results
final_region_table <- kw_results %>%
  filter(feature %in% sig_kw_features) %>%
  left_join(feature_pairwise_summary, by = "feature") %>%
  arrange(kw_p_adj_fdr)

# Regions significant across all 3 pairwise comparisons
regions_sig_all_3_pairs <- final_region_table %>%
  filter(sig_all_pairs) %>%
  pull(feature)

cat("Number of regions significant (KW FDR<0.05):", length(sig_kw_features), "\n")
cat("Number of regions significant in ALL 3 pairwise comparisons:", length(regions_sig_all_3_pairs), "\n")

# View top table
print(head(final_region_table, 30))

# If you want to save:
# write.csv(final_region_table, "KW_pairwise_region_results.csv", row.names = FALSE)


selected_features <- regions_sig_all_3_pairs
# or manually:
# selected_features <- c("Left.Hippocampus.Volume..ST17SV.", "Right.Entorhinal.Volume..ST83CV.", ...)
length(selected_features)

# Ensure DIAGNOSIS is a factor
data$DIAGNOSIS <- factor(data$DIAGNOSIS, levels = c("CN", "MCI", "AD"))

# Color mapping
diag_colors <- c(
  "CN"  = "darkgreen",
  "MCI" = "blue",
  "AD"  = "red"
)

# Numeric encoding for y-axis
y_pos <- as.numeric(data$DIAGNOSIS)

plot_feature_points <- function(feature_name) {
  
  x <- data[[feature_name]]
  
  # Remove NA safely
  keep <- !is.na(x)
  x <- x[keep]
  y <- y_pos[keep]
  col <- diag_colors[data$DIAGNOSIS[keep]]
  
  # Small vertical jitter
  y_jitter <- y + runif(length(y), -0.15, 0.15)
  
  plot(
    x, y_jitter,
    col = col,
    pch = 16,
    cex = 0.7,
    yaxt = "n",
    xlab = "Atrophy rate",
    ylab = "",
    main = feature_name
  )
  
  axis(2, at = 1:3, labels = c("CN", "MCI", "AD"), las = 1)
}
selected_features

par(mfrow = c(1, 1))  # ensure single plot at a time

plot_feature_points(selected_features[1])
plot_feature_points(selected_features[2])
plot_feature_points(selected_features[3])
plot_feature_points(selected_features[4])
plot_feature_points(selected_features[5])
plot_feature_points(selected_features[6])
plot_feature_points(selected_features[7])
plot_feature_points(selected_features[8])
plot_feature_points(selected_features[9])
plot_feature_points(selected_features[10])
plot_feature_points(selected_features[11])
plot_feature_points(selected_features[12])
plot_feature_points(selected_features[13])
plot_feature_points(selected_features[14])
plot_feature_points(selected_features[15])
plot_feature_points(selected_features[16])
plot_feature_points(selected_features[17])
plot_feature_points(selected_features[18])
plot_feature_points(selected_features[19])
plot_feature_points(selected_features[20])
plot_feature_points(selected_features[21])
plot_feature_points(selected_features[22])


## =========================
## INPUTS
## =========================

# Your selected features vector (length 22)
# selected_features <- regions_sig_all_3_pairs
# or whatever vector you already have

q_low  <- 0.05
q_high <- 0.95

# Make sure DIAGNOSIS is set and ordered
data$DIAGNOSIS <- factor(data$DIAGNOSIS, levels = c("CN", "MCI", "AD"))

# Colors
diag_colors <- c(CN = "darkgreen", MCI = "blue", AD = "red")

## =========================
## 1) Compute CN quantiles per feature
## =========================

cn_idx <- which(data$DIAGNOSIS == "CN")

# Matrices/vectors to store quantiles
Qlow  <- numeric(length(selected_features))
Qhigh <- numeric(length(selected_features))
names(Qlow) <- names(Qhigh) <- selected_features

for (j in seq_along(selected_features)) {
  feat <- selected_features[j]
  
  # If column missing, set NA and continue
  if (!feat %in% names(data)) {
    Qlow[j] <- NA
    Qhigh[j] <- NA
    next
  }
  
  x_cn <- data[cn_idx, feat]
  x_cn <- x_cn[!is.na(x_cn)]
  
  # If not enough data, set NA
  if (length(x_cn) < 5) {
    Qlow[j] <- NA
    Qhigh[j] <- NA
    next
  }
  
  qs <- quantile(x_cn, probs = c(q_low, q_high), na.rm = TRUE, type = 7)
  Qlow[j]  <- unname(qs[1])
  Qhigh[j] <- unname(qs[2])
}

# Optional: inspect quantile table
quantile_table <- data.frame(
  feature = selected_features,
  Qlow = Qlow,
  Qhigh = Qhigh
)
print(quantile_table)

## Remove features with missing thresholds (if any)
keep_feats <- which(!is.na(Qlow) & !is.na(Qhigh))
selected_features2 <- selected_features[keep_feats]
Qlow2  <- Qlow[keep_feats]
Qhigh2 <- Qhigh[keep_feats]

cat("Using", length(selected_features2), "features (after removing missing quantiles)\n")


## =========================
## 2) Compute violation matrix and per-patient violation counts
## =========================

n <- nrow(data)
p <- length(selected_features2)

viol_mat <- matrix(FALSE, nrow = n, ncol = p)
colnames(viol_mat) <- selected_features2

for (j in seq_along(selected_features2)) {
  feat <- selected_features2[j]
  x <- data[[feat]]
  
  # Violation if outside CN envelope; NA stays FALSE (or you can set NA)
  viol_mat[, j] <- (!is.na(x)) & (x < Qlow2[j] | x > Qhigh2[j])
}

# Count violations per patient
viol_count <- rowSums(viol_mat)

# Store back into data if you want
data$VIOL_COUNT <- viol_count


## =========================
## 3) Plot violation counts (points), like earlier
## =========================

# Jittered y positions for diagnosis
y_pos <- as.numeric(data$DIAGNOSIS)
y_jit <- y_pos + runif(n, -0.15, 0.15)
pt_col <- diag_colors[as.character(data$DIAGNOSIS)]

par(mfrow = c(1, 1))

plot(
  data$VIOL_COUNT, y_jit,
  col = pt_col,
  pch = 16, cex = 0.7,
  yaxt = "n",
  xlab = "Number of CN-quantile violations",
  ylab = "",
  main = paste0("Violations outside CN envelope (", q_low, "-", q_high, " quantiles)")
)
axis(2, at = 1:3, labels = c("CN", "MCI", "AD"), las = 1)


## =========================
## 4) OPTIONAL: one separate plot per class (CN / MCI / AD)
## =========================

par(mfrow = c(1, 1))

plot_one_class <- function(class_name) {
  idx <- which(data$DIAGNOSIS == class_name)
  x <- data$VIOL_COUNT[idx]
  
  # x-axis jitter so points don't pile up
  x_jit <- x + runif(length(x), -0.15, 0.15)
  
  plot(
    x_jit, rep(1, length(x_jit)),
    col = diag_colors[class_name],
    pch = 16, cex = 0.8,
    yaxt = "n", ylab = "",
    xlab = "Number of CN-quantile violations",
    main = paste0(class_name, " patients: violation counts")
  )
}

plot_one_class("CN")
plot_one_class("MCI")
plot_one_class("AD")

# Diagnosis levels (ensure order)
classes <- c("CN", "MCI", "AD")

# Initialize result matrix
viol_summary <- matrix(
  0,
  nrow = length(classes),
  ncol = length(selected_features2)
)

rownames(viol_summary) <- classes
colnames(viol_summary) <- selected_features2

# Fill the table
for (i in seq_along(classes)) {
  class_idx <- which(data$DIAGNOSIS == classes[i])
  
  # Count violations per feature within class
  viol_summary[i, ] <- colSums(viol_mat[class_idx, , drop = FALSE])
}

# Convert to data frame for easy viewing
viol_summary_df <- as.data.frame(viol_summary)

viol_summary_df

#####

# MCI patients sorted by number of violations
mci_sorted <- data[data$DIAGNOSIS == "MCI", c("RID", "GENOTYPE", "VIOL_COUNT")]
mci_sorted <- mci_sorted[order(mci_sorted$VIOL_COUNT, decreasing = TRUE), ]
row.names(mci_sorted) <- NULL

# AD patients sorted by number of violations
ad_sorted <- data[data$DIAGNOSIS == "AD", c("RID", "GENOTYPE", "VIOL_COUNT")]
ad_sorted <- ad_sorted[order(ad_sorted$VIOL_COUNT, decreasing = TRUE), ]
row.names(ad_sorted) <- NULL

# View them
head(mci_sorted, 20)
head(ad_sorted, 20)

tail(mci_sorted, 80)
tail(ad_sorted, 40)

#################

## Extract violation counts
mci_counts <- data$VIOL_COUNT[data$DIAGNOSIS == "MCI"]
ad_counts  <- data$VIOL_COUNT[data$DIAGNOSIS == "AD"]

## Create y positions
y_mci <- rep(1, length(mci_counts))
y_ad  <- rep(2, length(ad_counts))

## Jitter for readability
y_mci_j <- y_mci + runif(length(y_mci), -0.12, 0.12)
y_ad_j  <- y_ad  + runif(length(y_ad),  -0.12, 0.12)

## Plot
plot(
  mci_counts, y_mci_j,
  col = "blue",
  pch = 16, cex = 0.7,
  xlim = range(c(mci_counts, ad_counts)),
  ylim = c(0.5, 2.5),
  yaxt = "n",
  xlab = "Number of CN-quantile violations",
  ylab = "",
  main = "Violation counts per patient (MCI vs AD)"
)

points(
  ad_counts, y_ad_j,
  col = "red",
  pch = 16, cex = 0.7
)

axis(2, at = c(1, 2), labels = c("MCI", "AD"), las = 1)

legend(
  "topright",
  legend = c("MCI", "AD"),
  col = c("blue", "red"),
  pch = 16,
  bty = "n"
)

#####################àà

## Extract violation counts
mci_counts <- data$VIOL_COUNT[data$DIAGNOSIS == "MCI"]
ad_counts  <- data$VIOL_COUNT[data$DIAGNOSIS == "AD"]

## Create x positions (just to spread points)
x_mci <- seq_along(mci_counts)
x_ad  <- seq_along(ad_counts)

## Combine limits
y_lim <- range(c(mci_counts, ad_counts))

## Plot MCI first
plot(
  x_mci, mci_counts,
  col = "blue",
  pch = 16,
  cex = 0.7,
  ylim = y_lim,
  xlab = "Patient index",
  ylab = "Total number of CN-quantile violations",
  main = "Violation count per patient (MCI = blue, AD = red)"
)

## Add AD points
points(
  x_ad, ad_counts,
  col = "red",
  pch = 16,
  cex = 0.7
)

legend(
  "topleft",
  legend = c("MCI", "AD"),
  col = c("blue", "red"),
  pch = 16,
  bty = "n"
)

#### BOXPLOT

# Extract violation counts
mci_counts <- data$VIOL_COUNT[data$DIAGNOSIS == "MCI"]
ad_counts  <- data$VIOL_COUNT[data$DIAGNOSIS == "AD"]

# Combine into a list for boxplot
viol_list <- list(
  MCI = mci_counts,
  AD  = ad_counts
)

# Draw boxplot
boxplot(
  viol_list,
  col = c("blue", "red"),
  border = c("blue", "red"),
  ylab = "Total number of CN-quantile violations",
  main = "Violation count per patient (MCI vs AD)",
  outline = TRUE
)

# Optional: add grid for readability
grid(nx = NA, ny = NULL)


######### CLASSIFIER

# Initialize predicted diagnosis
pred_diag <- character(nrow(data))

# Apply the rule
pred_diag[data$VIOL_COUNT <= 2] <- "CN"
pred_diag[data$VIOL_COUNT >= 3 & data$VIOL_COUNT <= 9] <- "MCI"
pred_diag[data$VIOL_COUNT > 9] <- "AD"

# Convert to factor with same levels
pred_diag <- factor(pred_diag, levels = c("CN", "MCI", "AD"))

# Store if you want
data$PRED_DIAG <- pred_diag


conf_mat <- table(
  True = data$DIAGNOSIS,
  Predicted = data$PRED_DIAG
)

conf_mat

overall_accuracy <- sum(diag(conf_mat)) / sum(conf_mat)
cat("Overall accuracy:", round(overall_accuracy, 3), "\n")

class_accuracy <- diag(conf_mat) / rowSums(conf_mat)
round(class_accuracy, 3)

balanced_accuracy <- mean(class_accuracy, na.rm = TRUE)
cat("Balanced accuracy:", round(balanced_accuracy, 3), "\n")

################ BINARY CLASSIFIER

# Choose threshold
t_AD <- 6  # <- change this to try 6,7,8,10,...

# Initialize predicted label
pred_bin <- character(nrow(data))

# Apply rule
pred_bin[data$VIOL_COUNT >= t_AD] <- "AD"
pred_bin[data$VIOL_COUNT <  t_AD] <- "non-AD"

# Convert to factor
pred_bin <- factor(pred_bin, levels = c("non-AD", "AD"))

# True binary labels
true_bin <- ifelse(data$DIAGNOSIS == "AD", "AD", "non-AD")
true_bin <- factor(true_bin, levels = c("non-AD", "AD"))

conf_mat_bin <- table(
  True = true_bin,
  Predicted = pred_bin
)

conf_mat_bin

overall_accuracy <- sum(diag(conf_mat_bin)) / sum(conf_mat_bin)
cat("Overall accuracy:", round(overall_accuracy, 3), "\n")

sensitivity_AD <- conf_mat_bin["AD","AD"] /
  sum(conf_mat_bin["AD",])

specificity_nonAD <- conf_mat_bin["non-AD","non-AD"] /
  sum(conf_mat_bin["non-AD",])

cat("Specificity (non-AD recall):", round(specificity_nonAD, 3), "\n")
cat("Sensitivity (AD recall):", round(sensitivity_AD, 3), "\n")


balanced_accuracy <- (sensitivity_AD + specificity_nonAD) / 2
cat("Balanced accuracy:", round(balanced_accuracy, 3), "\n")

# ---- BINARY CLASSIFIER WITH MAGNITUDE ----

## =========================================================
## Nested envelope scoring: CN envelope (0), outside CN (1),
## outside MCI (2). Then classify AD vs non-AD by threshold T.
## =========================================================

# PARAMETERS
q_low  <- 0.05
q_high <- 0.95

T <- 9   # <-- change this threshold and re-run

# Ensure diagnosis factor
data$DIAGNOSIS <- factor(data$DIAGNOSIS, levels = c("CN", "MCI", "AD"))

# Indices
idx_cn  <- which(data$DIAGNOSIS == "CN")
idx_mci <- which(data$DIAGNOSIS == "MCI")

p <- length(selected_features2)
n <- nrow(data)

# Storage for envelopes
Qlow_CN  <- numeric(p); Qhigh_CN  <- numeric(p)
Qlow_MCI <- numeric(p); Qhigh_MCI <- numeric(p)

names(Qlow_CN) <- names(Qhigh_CN) <- selected_features2
names(Qlow_MCI) <- names(Qhigh_MCI) <- selected_features2

# Compute CN and MCI quantiles per feature
for (j in seq_along(selected_features2)) {
  feat <- selected_features2[j]
  if (!feat %in% names(data)) {
    Qlow_CN[j] <- NA; Qhigh_CN[j] <- NA
    Qlow_MCI[j] <- NA; Qhigh_MCI[j] <- NA
    next
  }
  
  x_cn <- data[idx_cn, feat]
  x_cn <- x_cn[!is.na(x_cn)]
  x_mci <- data[idx_mci, feat]
  x_mci <- x_mci[!is.na(x_mci)]
  
  if (length(x_cn) < 5 || length(x_mci) < 5) {
    Qlow_CN[j] <- NA; Qhigh_CN[j] <- NA
    Qlow_MCI[j] <- NA; Qhigh_MCI[j] <- NA
    next
  }
  
  qs_cn  <- quantile(x_cn,  probs = c(q_low, q_high), na.rm = TRUE, type = 7)
  qs_mci <- quantile(x_mci, probs = c(q_low, q_high), na.rm = TRUE, type = 7)
  
  Qlow_CN[j]  <- unname(qs_cn[1]);  Qhigh_CN[j]  <- unname(qs_cn[2])
  Qlow_MCI[j] <- unname(qs_mci[1]); Qhigh_MCI[j] <- unname(qs_mci[2])
}

# Drop features with missing envelopes
keep <- which(!is.na(Qlow_CN) & !is.na(Qhigh_CN) & !is.na(Qlow_MCI) & !is.na(Qhigh_MCI))
feats <- selected_features2[keep]
Qlow_CN  <- Qlow_CN[keep];  Qhigh_CN  <- Qhigh_CN[keep]
Qlow_MCI <- Qlow_MCI[keep]; Qhigh_MCI <- Qhigh_MCI[keep]

cat("Using", length(feats), "features (after dropping missing envelopes)\n")

# Build nested score matrix: 0/1/2
score_mat <- matrix(0, nrow = n, ncol = length(feats))
colnames(score_mat) <- feats

for (j in seq_along(feats)) {
  feat <- feats[j]
  x <- data[[feat]]
  
  ok <- !is.na(x)
  
  out_CN  <- ok & (x < Qlow_CN[j]  | x > Qhigh_CN[j])
  out_MCI <- ok & (x < Qlow_MCI[j] | x > Qhigh_MCI[j])
  
  # default 0
  score_mat[out_CN,  j] <- 1
  score_mat[out_MCI, j] <- 2
}

# Total score per patient
data$NESTED_SCORE <- rowSums(score_mat)

## -------------------------
## Binary classification: AD vs non-AD
## -------------------------
true_bin <- ifelse(data$DIAGNOSIS == "AD", "AD", "non-AD")
true_bin <- factor(true_bin, levels = c("non-AD", "AD"))

pred_bin <- ifelse(data$NESTED_SCORE >= T, "AD", "non-AD")
pred_bin <- factor(pred_bin, levels = c("non-AD", "AD"))

conf_mat_bin <- table(True = true_bin, Predicted = pred_bin)
conf_mat_bin

overall_accuracy <- sum(diag(conf_mat_bin)) / sum(conf_mat_bin)

sensitivity_AD <- conf_mat_bin["AD","AD"] / sum(conf_mat_bin["AD",])
specificity_nonAD <- conf_mat_bin["non-AD","non-AD"] / sum(conf_mat_bin["non-AD",])

balanced_accuracy <- (sensitivity_AD + specificity_nonAD) / 2

cat("q_low =", q_low, " q_high =", q_high, " Threshold T =", T, "\n")
cat("Overall accuracy:", round(overall_accuracy, 3), "\n")
cat("Specificity (non-AD recall):", round(specificity_nonAD, 3), "\n")
cat("Sensitivity (AD recall):", round(sensitivity_AD, 3), "\n")
cat("Balanced accuracy:", round(balanced_accuracy, 3), "\n")


######## FULL CLASSIFIER

## =========================================================
## STAGE 1: AD vs non-AD (already validated)
## =========================================================

t_AD <- 7   # fixed from your previous results

pred_stage1 <- character(nrow(data))
pred_stage1[data$VIOL_COUNT >= t_AD] <- "AD"
pred_stage1[data$VIOL_COUNT <  t_AD] <- "non-AD"

pred_stage1 <- factor(pred_stage1, levels = c("non-AD", "AD"))

data$PRED_STAGE1 <- pred_stage1


## =========================================================
## STAGE 2: CN vs MCI (only among predicted non-AD)
## =========================================================

# Parameters for CN envelope
q_low2  <- 0.05
q_high2 <- 0.95

# Threshold for MCI classification
t_MCI <- 3   # <- YOU WILL TUNE THIS (try 1,2,3,...)

# Subset: only predicted non-AD
idx_nonAD <- which(data$PRED_STAGE1 == "non-AD")
data_nonAD <- data[idx_nonAD, ]

# True labels here are CN / MCI only
table(data_nonAD$DIAGNOSIS)

# Indices of CN patients (within non-AD)
idx_cn2 <- which(data_nonAD$DIAGNOSIS == "CN")

# Compute CN quantiles per feature
p <- length(selected_features2)

Qlow_CN2  <- numeric(p)
Qhigh_CN2 <- numeric(p)

names(Qlow_CN2) <- names(Qhigh_CN2) <- selected_features2

for (j in seq_along(selected_features2)) {
  feat <- selected_features2[j]
  
  x_cn <- data_nonAD[idx_cn2, feat]
  x_cn <- x_cn[!is.na(x_cn)]
  
  if (length(x_cn) < 5) {
    Qlow_CN2[j] <- NA
    Qhigh_CN2[j] <- NA
    next
  }
  
  qs <- quantile(x_cn, probs = c(q_low2, q_high2), na.rm = TRUE)
  Qlow_CN2[j]  <- unname(qs[1])
  Qhigh_CN2[j] <- unname(qs[2])
}

# Keep usable features
keep2 <- which(!is.na(Qlow_CN2) & !is.na(Qhigh_CN2))
feats2 <- selected_features2[keep2]
Qlow_CN2  <- Qlow_CN2[keep2]
Qhigh_CN2 <- Qhigh_CN2[keep2]

cat("Stage 2 uses", length(feats2), "features\n")

# Count CN-envelope violations
viol2_mat <- matrix(0, nrow = nrow(data_nonAD), ncol = length(feats2))
colnames(viol2_mat) <- feats2

for (j in seq_along(feats2)) {
  x <- data_nonAD[[feats2[j]]]
  viol2_mat[, j] <- (!is.na(x)) & (x < Qlow_CN2[j] | x > Qhigh_CN2[j])
}

data_nonAD$VIOL_COUNT_STAGE2 <- rowSums(viol2_mat)

# Stage-2 prediction
pred_stage2 <- character(nrow(data_nonAD))
pred_stage2[data_nonAD$VIOL_COUNT_STAGE2 >= t_MCI] <- "MCI"
pred_stage2[data_nonAD$VIOL_COUNT_STAGE2 <  t_MCI] <- "CN"

pred_stage2 <- factor(pred_stage2, levels = c("CN", "MCI"))

data_nonAD$PRED_STAGE2 <- pred_stage2


## =========================================================
## COMBINE STAGE 1 + STAGE 2
## =========================================================

final_pred <- character(nrow(data))
final_pred[data$PRED_STAGE1 == "AD"] <- "AD"
final_pred[idx_nonAD] <- as.character(data_nonAD$PRED_STAGE2)

final_pred <- factor(final_pred, levels = c("CN", "MCI", "AD"))
data$PRED_FINAL <- final_pred


## =========================================================
## EVALUATION: 3-CLASS PERFORMANCE
## =========================================================

# Confusion matrix
conf_mat_3 <- table(
  True = data$DIAGNOSIS,
  Predicted = data$PRED_FINAL
)

conf_mat_3

# Overall accuracy
overall_accuracy <- sum(diag(conf_mat_3)) / sum(conf_mat_3)
cat("Overall accuracy:", round(overall_accuracy, 3), "\n")

# Per-class recall
class_recall <- diag(conf_mat_3) / rowSums(conf_mat_3)
round(class_recall, 3)

# Balanced accuracy (mean recall)
balanced_accuracy <- mean(class_recall, na.rm = TRUE)
cat("Balanced accuracy:", round(balanced_accuracy, 3), "\n")















