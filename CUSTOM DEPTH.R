#---- 0) Load cleaned data ----
data <- read.csv("dataQR_clean.csv")

# Identify feature columns (exclude ID + label columns)
feature_cols <- setdiff(
  colnames(data),
  c("RID", "DIAGNOSIS", "GENOTYPE", "DIAGNOSIS_NUM")
)

# Ensure numeric
data[feature_cols] <- lapply(data[feature_cols], as.numeric)

# Robust center and scale
medians <- apply(data[feature_cols], 2, median, na.rm = TRUE)
mads    <- apply(data[feature_cols], 2, mad,    na.rm = TRUE)

# Avoid division by zero
mads[mads == 0] <- 1e-6

# Compute robust z-scores
robust_z <- sweep(data[feature_cols], 2, medians, "-")
robust_z <- sweep(robust_z, 2, mads, "/")

# Threshold for extreme outliers
threshold <- 188885   # 4 MAD is VERY extreme; can increase to 10 if needed

# Identify rows with any extreme feature
extreme_rows <- apply(abs(robust_z) > threshold, 1, any)

cat("Number of rows removed:", sum(extreme_rows), "\n")

# Remove those rows
data_clean <- data[!extreme_rows, ]

cat("Remaining rows:", nrow(data_clean), "\n")

## =========================================================
## CN-baseline nested envelopes + two-term "depth" components
##  Term 1: ordinal envelope violations (0 / 1 / 3)
##  Term 2: global CN-deviation energy (sum of squared z)
##  Outputs: per-class mean + SD for each term
## =========================================================

library(dplyr)

#----------------------------
# 0) Inputs you already have
#----------------------------
# data_clean : your cleaned dataset (rows with extreme outliers removed)
# feature_cols : vector of feature column names
# data_clean$DIAGNOSIS : factor with levels c("CN","MCI","AD") (order not critical)

# If you haven't defined feature_cols yet, uncomment:
# feature_cols <- setdiff(colnames(data_clean), c("RID","DIAGNOSIS","GENOTYPE","DIAGNOSIS_NUM"))

# Ensure numeric
data_clean[feature_cols] <- lapply(data_clean[feature_cols], as.numeric)

#----------------------------
# 1) Define CN baseline + nested envelopes (CN inside MCI)
#----------------------------
# Using CN median + CN MAD (robust). Envelopes are:
#   CN band : m ± a * s
#   MCI band: m ± b * s , with b > a
a <- 1.0   # tighter band (CN)
b <- 2.5   # wider band (MCI); increase if you want fewer "AD/outside" flags

# Ordinal penalties (tune if needed)
pen_cn  <- 0
pen_mci <- 1
pen_ad  <- 3

cn_df <- data_clean %>% filter(DIAGNOSIS == "CN")

m <- apply(cn_df[, feature_cols, drop = FALSE], 2, median, na.rm = TRUE)
s <- apply(cn_df[, feature_cols, drop = FALSE], 2, mad,    na.rm = TRUE)
s[s == 0] <- 1e-6

cn_low  <- m - a * s
cn_high <- m + a * s
mci_low  <- m - b * s
mci_high <- m + b * s

#----------------------------
# 2) Compute Term 1 and Term 2 for each subject
#----------------------------
X <- as.matrix(data_clean[, feature_cols, drop = FALSE])

# Term 1: per-feature category relative to nested envelopes
inside_cn  <- sweep(X, 2, cn_low,  `>=`) & sweep(X, 2, cn_high, `<=`)
inside_mci <- sweep(X, 2, mci_low, `>=`) & sweep(X, 2, mci_high, `<=`)

# Per-feature penalty: CN=0, MCI-ring=1, Outside=3
# (MCI-ring = inside_mci but not inside_cn)
penalty_mat <- matrix(pen_ad, nrow = nrow(X), ncol = ncol(X))
penalty_mat[inside_mci] <- pen_mci
penalty_mat[inside_cn]  <- pen_cn

term1_raw <- rowSums(penalty_mat, na.rm = TRUE)

# Optional: normalize by number of features (recommended so scale doesn't depend on p)
p <- length(feature_cols)
term1 <- term1_raw / p

# Term 2: CN-centered squared robust z-score sum (energy)
Z <- sweep(X, 2, m, "-")
Z <- sweep(Z, 2, s, "/")
term2_raw <- rowSums(Z^2, na.rm = TRUE)

# Optional: normalize by p (recommended)
term2 <- term2_raw / p

# Attach to dataset
scored <- data_clean %>%
  mutate(
    term1 = term1,
    term2 = term2
  )

#----------------------------
# 3) Per-class summaries (mean + SD)
#----------------------------

# Ensure desired order
scored$DIAGNOSIS <- factor(
  scored$DIAGNOSIS,
  levels = c("CN", "MCI", "AD")
)

summary_terms <- scored %>%
  group_by(DIAGNOSIS) %>%
  summarise(
    n = n(),
    term1_mean = mean(term1, na.rm = TRUE),
    term1_sd   = sd(term1,   na.rm = TRUE),
    term2_mean = mean(term2, na.rm = TRUE),
    term2_sd   = sd(term2,   na.rm = TRUE),
    .groups = "drop"
  )

print(summary_terms)

library(ggplot2)

# Ensure class order
scored$DIAGNOSIS <- factor(scored$DIAGNOSIS, levels = c("CN","MCI","AD"))

# --- Boxplot for Term 1 ---
p1 <- ggplot(scored, aes(x = DIAGNOSIS, y = term1)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5) +
  labs(
    title = "Term 1 (Envelope violation score) by diagnosis",
    x = "Diagnosis",
    y = "Term 1 (normalized)"
  ) +
  theme_bw(base_size = 14)

print(p1)

# --- Boxplot for Term 2 ---
p2 <- ggplot(scored, aes(x = DIAGNOSIS, y = term2)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5) +
  labs(
    title = "Term 2 (CN-deviation energy) by diagnosis",
    x = "Diagnosis",
    y = "Term 2 (normalized)"
  ) +
  theme_bw(base_size = 14)

print(p2)

p2_log <- ggplot(scored, aes(x = DIAGNOSIS, y = term2)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5) +
  scale_y_log10() +
  labs(
    title = "Term 2 (CN-deviation energy) by diagnosis (log scale)",
    x = "Diagnosis",
    y = "Term 2 (normalized, log10 scale)"
  ) +
  theme_bw(base_size = 14)

print(p2_log)

## Outliergram

library(ggplot2)
library(patchwork)

# Ensure correct order
scored$DIAGNOSIS <- factor(scored$DIAGNOSIS, levels = c("CN","MCI","AD"))

# Color mapping
cols <- c("CN" = "green", "MCI" = "blue", "AD" = "red")

# --- Individual class plots ---
p_cn <- ggplot(subset(scored, DIAGNOSIS == "CN"),
               aes(x = term1, y = term2)) +
  geom_point(color = "green", alpha = 0.6) +
  labs(title = "CN", x = "Term 1", y = "Term 2") +
  theme_bw(base_size = 14)

p_mci <- ggplot(subset(scored, DIAGNOSIS == "MCI"),
                aes(x = term1, y = term2)) +
  geom_point(color = "blue", alpha = 0.6) +
  labs(title = "MCI", x = "Term 1", y = "Term 2") +
  theme_bw(base_size = 14)

p_ad <- ggplot(subset(scored, DIAGNOSIS == "AD"),
               aes(x = term1, y = term2)) +
  geom_point(color = "red", alpha = 0.6) +
  labs(title = "AD", x = "Term 1", y = "Term 2") +
  theme_bw(base_size = 14)

# --- Combined plot ---
p_all <- ggplot(scored,
                aes(x = term1, y = term2, color = DIAGNOSIS)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = cols) +
  labs(title = "All Classes",
       x = "Term 1",
       y = "Term 2",
       color = "Diagnosis") +
  theme_bw(base_size = 14)

# --- Arrange in 2x2 layout ---
(p_cn | p_mci) /
  (p_ad | p_all)



library(ggplot2)
library(patchwork)

# Ensure correct order
scored$DIAGNOSIS <- factor(scored$DIAGNOSIS, levels = c("CN","MCI","AD"))

# Color mapping
cols <- c("CN" = "green", "MCI" = "blue", "AD" = "red")

# Ensure term2 is strictly positive for log scale
if (any(scored$term2 <= 0)) {
  shift_val <- min(scored$term2[scored$term2 > 0]) * 0.5
  scored$term2_plot <- scored$term2 + shift_val
} else {
  scored$term2_plot <- scored$term2
}

# --- Individual class plots ---
p_cn <- ggplot(subset(scored, DIAGNOSIS == "CN"),
               aes(x = term1, y = term2_plot)) +
  geom_point(color = "green", alpha = 0.6) +
  scale_y_log10() +
  labs(title = "CN", x = "Term 1", y = "Term 2 (log scale)") +
  theme_bw(base_size = 14)

p_mci <- ggplot(subset(scored, DIAGNOSIS == "MCI"),
                aes(x = term1, y = term2_plot)) +
  geom_point(color = "blue", alpha = 0.6) +
  scale_y_log10() +
  labs(title = "MCI", x = "Term 1", y = "Term 2 (log scale)") +
  theme_bw(base_size = 14)

p_ad <- ggplot(subset(scored, DIAGNOSIS == "AD"),
               aes(x = term1, y = term2_plot)) +
  geom_point(color = "red", alpha = 0.6) +
  scale_y_log10() +
  labs(title = "AD", x = "Term 1", y = "Term 2 (log scale)") +
  theme_bw(base_size = 14)

# --- Combined plot ---
p_all <- ggplot(scored,
                aes(x = term1, y = term2_plot, color = DIAGNOSIS)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = cols) +
  scale_y_log10() +
  labs(title = "All Classes",
       x = "Term 1",
       y = "Term 2 (log scale)",
       color = "Diagnosis") +
  theme_bw(base_size = 14)

# --- Arrange in 2x2 layout ---
(p_cn | p_mci) /
  (p_ad | p_all)


# CORRELATION

cor_total <- cor(scored$term1, scored$term2, use = "complete.obs")
cat("Overall correlation (Term1 vs Term2):", cor_total, "\n")

library(dplyr)

cor_by_class <- scored %>%
  group_by(DIAGNOSIS) %>%
  summarise(
    correlation = cor(term1, term2, use = "complete.obs"),
    .groups = "drop"
  )

print(cor_by_class)


## BOXPLOT OF THE SUM

library(ggplot2)

# Ensure ordering
scored$DIAGNOSIS <- factor(scored$DIAGNOSIS, levels = c("CN","MCI","AD"))

# Color mapping
cols <- c("CN" = "green", "MCI" = "blue", "AD" = "red")

# Combined score 1
scored$sum_raw <- scored$term1 + scored$term2

p_sum_raw <- ggplot(scored, aes(x = DIAGNOSIS, y = sum_raw, fill = DIAGNOSIS)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = cols) +
  labs(
    title = "Boxplot of Term1 + Term2",
    x = "Diagnosis",
    y = "Term1 + Term2"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

print(p_sum_raw)

## BOXPLOT OF THE LOG SUM

# Ensure positive for log
if (any(scored$term2 <= 0)) {
  shift_val <- min(scored$term2[scored$term2 > 0]) * 0.5
  scored$term2_log <- log(scored$term2 + shift_val)
} else {
  scored$term2_log <- log(scored$term2)
}

# Combined score 2
scored$sum_log <- scored$term1 + scored$term2_log

p_sum_log <- ggplot(scored, aes(x = DIAGNOSIS, y = sum_log, fill = DIAGNOSIS)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = cols) +
  labs(
    title = "Boxplot of Term1 + log(Term2)",
    x = "Diagnosis",
    y = "Term1 + log(Term2)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

print(p_sum_log)


## TERM 1 TUNING

## =========================================================
## Grid search for envelope penalties with:
##  (1) progress bar
##  (2) only try w2 > w1
##  (3) boxplot for the best combination
## =========================================================

library(dplyr)
library(ggplot2)
library(progress)

# -----------------------------
# Assumes you already have:
#   X           : matrix of features (n x p)
#   p           : number of features (length(feature_cols))
#   inside_cn   : logical matrix (n x p)
#   inside_mci  : logical matrix (n x p)
#   scored      : data.frame with DIAGNOSIS and term2_log
# -----------------------------

# Make sure order is consistent
scored$DIAGNOSIS <- factor(scored$DIAGNOSIS, levels = c("CN","MCI","AD"))
cols <- c("CN"="green", "MCI"="blue", "AD"="red")

# ---- Large grids (edit freely) ----
# w1 = penalty for MCI ring, w2 = penalty for outside (AD region)
w1_vals <- seq(0, 30, by = 0.1)
w2_vals <- seq(0, 30, by = 0.1)

# ---- Build only pairs with w2 > w1 ----
grid <- expand.grid(w1 = w1_vals, w2 = w2_vals) %>%
  filter(w2 > w1)

cat("Number of (w1,w2) pairs to evaluate:", nrow(grid), "\n")

# ---- Progress bar ----
pb <- progress_bar$new(
  format = "Grid search [:bar] :current/:total (:percent) | ETA: :eta",
  total  = nrow(grid),
  clear  = FALSE,
  width  = 70
)

# ---- Separation metric storage ----
sep <- numeric(nrow(grid))

# ---- Loop ----
for (i in seq_len(nrow(grid))) {
  w1 <- grid$w1[i]
  w2 <- grid$w2[i]
  
  # Term 1 with penalties:
  #   inside CN -> 0
  #   inside MCI ring (inside_mci & !inside_cn) -> w1
  #   outside MCI -> w2
  penalty_mat <- matrix(w2, nrow = nrow(X), ncol = ncol(X))
  penalty_mat[inside_mci] <- w1
  penalty_mat[inside_cn]  <- 0
  
  term1_new <- rowSums(penalty_mat) / p
  
  # Combined score: Term1 + log(Term2)
  score <- term1_new + scored$term2_log
  
  # Between/within separation ratio
  group_means  <- tapply(score, scored$DIAGNOSIS, mean, na.rm = TRUE)
  overall_mean <- mean(score, na.rm = TRUE)
  
  n_g <- table(scored$DIAGNOSIS)
  
  between_var <- sum(n_g * (group_means - overall_mean)^2, na.rm = TRUE)
  within_var  <- sum(tapply(score, scored$DIAGNOSIS, var, na.rm = TRUE), na.rm = TRUE)
  
  sep[i] <- between_var / within_var
  
  pb$tick()
}

grid$separation <- sep

# ---- Best combination ----
best_idx <- which.max(grid$separation)
best <- grid[best_idx, ]
print(best)

best_w1 <- best$w1
best_w2 <- best$w2

# ---- Recompute best score and plot boxplot ----
penalty_mat_best <- matrix(best_w2, nrow = nrow(X), ncol = ncol(X))
penalty_mat_best[inside_mci] <- best_w1
penalty_mat_best[inside_cn]  <- 0

term1_best <- rowSums(penalty_mat_best) / p
score_best <- term1_best + scored$term2_log

plot_df <- scored %>%
  mutate(score_best = score_best)

p_best <- ggplot(plot_df, aes(x = DIAGNOSIS, y = score_best, fill = DIAGNOSIS)) +
  geom_boxplot(alpha = 0.75, outlier.size = 1.2) +
  scale_fill_manual(values = cols) +
  labs(
    title = paste0("Best penalties: w1 = ", best_w1, ", w2 = ", best_w2),
    subtitle = paste0("Separation ratio = ", round(best$separation, 4)),
    x = "Diagnosis",
    y = "Term1(w1,w2) + log(Term2)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

print(p_best)

# ---- (Optional) keep best score + best term1 in your main data ----
scored$term1_best <- term1_best
scored$score_best <- score_best

cor_overall <- cor(scored$term1_best,
                   scored$term2_log,
                   use = "complete.obs")

cat("Overall correlation (term1_best vs log_term2):",
    round(cor_overall, 4), "\n")

library(ggplot2)

ggplot(scored, aes(x = term1_best,
                   y = term2_log,
                   color = DIAGNOSIS)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("CN"="green",
                                "MCI"="blue",
                                "AD"="red")) +
  theme_bw(base_size = 14) +
  labs(title = "Correlation: Term1_best vs log(Term2)",
       x = "Term1 (best penalties)",
       y = "log(Term2)")




## =========================================================
## Option 1 + 2: build alternative Term2s and check correlation
##   Option 1: L1 deviation from CN baseline  -> term2_L1_log
##   Option 2: variance profile of deviations -> term3_var_absZ
## Compare both to term1_best
## =========================================================

library(dplyr)

# Assumes you already have:
#  - scored (data.frame) with: DIAGNOSIS, term1_best
#  - feature_cols (character vector)
#  - data_clean (or the same rows used to build scored)
#  - CN baseline center m and scale s (or we recompute them here)

# ----------------------------
# 0) Recompute CN baseline (safe)
# ----------------------------
scored$DIAGNOSIS <- factor(scored$DIAGNOSIS, levels = c("CN","MCI","AD"))

# Make sure feature matrix aligns with scored rows
# (If scored was created from data_clean in the same order, this is fine.)
X <- as.matrix(data_clean[, feature_cols, drop = FALSE])

cn_df <- data_clean[data_clean$DIAGNOSIS == "CN", , drop = FALSE]
m <- apply(cn_df[, feature_cols, drop = FALSE], 2, median, na.rm = TRUE)
s <- apply(cn_df[, feature_cols, drop = FALSE], 2, mad,    na.rm = TRUE)
s[s == 0] <- 1e-6

# Robust z-scores vs CN baseline
Z <- sweep(X, 2, m, "-")
Z <- sweep(Z, 2, s, "/")

p <- length(feature_cols)

# ----------------------------
# 1) Option 1: L1 deviation (then log)
# ----------------------------
term2_L1 <- rowSums(abs(Z), na.rm = TRUE) / p

# log needs positive values; add tiny epsilon just in case
eps <- 1e-12
term2_L1_log <- log(term2_L1 + eps)

# ----------------------------
# 2) Option 2: Variance profile of deviations
#    (captures "how concentrated vs spread" the deviations are)
# ----------------------------
term3_var_absZ <- apply(abs(Z), 1, var, na.rm = TRUE)

# Attach to scored
scored$term2_L1_log    <- term2_L1_log
scored$term3_var_absZ  <- term3_var_absZ

# ----------------------------
# 3) Correlations with term1_best
# ----------------------------
cor_term1_L1  <- cor(scored$term1_best, scored$term2_L1_log,   use = "complete.obs")
cor_term1_var <- cor(scored$term1_best, scored$term3_var_absZ, use = "complete.obs")

cat("Correlation(term1_best, term2_L1_log):   ", round(cor_term1_L1, 4),  "\n")
cat("Correlation(term1_best, term3_var_absZ): ", round(cor_term1_var, 4), "\n")

# ----------------------------
# 4) (Optional) correlations by class
# ----------------------------
cor_by_class <- scored %>%
  group_by(DIAGNOSIS) %>%
  summarise(
    cor_term1_L1  = cor(term1_best, term2_L1_log,   use = "complete.obs"),
    cor_term1_var = cor(term1_best, term3_var_absZ, use = "complete.obs"),
    .groups = "drop"
  )

print(cor_by_class)

# ----------------------------
# 5) (Optional) full correlation matrix
# ----------------------------
cor_mat <- cor(
  scored[, c("term1_best", "term2_L1_log", "term3_var_absZ")],
  use = "complete.obs"
)
print(round(cor_mat, 4))



### PLOT TERM 3

library(ggplot2)

# Ensure class order
scored$DIAGNOSIS <- factor(scored$DIAGNOSIS, levels = c("CN","MCI","AD"))

# Colors
cols <- c("CN"="green", "MCI"="blue", "AD"="red")

p_term3 <- ggplot(scored, aes(x = DIAGNOSIS, y = term3_var_absZ, fill = DIAGNOSIS)) +
  geom_boxplot(alpha = 0.75, outlier.size = 1.2) +
  scale_fill_manual(values = cols) +
  labs(
    title = "Variance of Absolute CN-Deviations by Diagnosis",
    x = "Diagnosis",
    y = "Variance of |Z|"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

print(p_term3)

# LOG VERSION TERM 3

p_term3_log <- ggplot(scored, aes(x = DIAGNOSIS, y = term3_var_absZ, fill = DIAGNOSIS)) +
  geom_boxplot(alpha = 0.75, outlier.size = 1.2) +
  scale_fill_manual(values = cols) +
  scale_y_log10() +
  labs(
    title = "Variance of |Z| by Diagnosis (log scale)",
    x = "Diagnosis",
    y = "Variance of |Z| (log scale)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

print(p_term3_log)

# Ensure L1 term exists
scored$term1_L1 <- term2_L1  # from previous computation (mean |Z|)

library(ggplot2)

scored$DIAGNOSIS <- factor(scored$DIAGNOSIS, levels = c("CN","MCI","AD"))
cols <- c("CN"="green", "MCI"="blue", "AD"="red")

p_L1 <- ggplot(scored, aes(x = DIAGNOSIS, y = term1_L1, fill = DIAGNOSIS)) +
  geom_boxplot(alpha = 0.75, outlier.size = 1.2) +
  scale_fill_manual(values = cols) +
  labs(
    title = "Mean Absolute CN-Deviation (L1) by Diagnosis",
    x = "Diagnosis",
    y = "Mean |Z|"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

print(p_L1)


# Combined score
scored$score_L1_var <- scored$term1_L1 + scored$term3_var_absZ

p_combined <- ggplot(scored, aes(x = DIAGNOSIS, y = score_L1_var, fill = DIAGNOSIS)) +
  geom_boxplot(alpha = 0.75, outlier.size = 1.2) +
  scale_fill_manual(values = cols) +
  labs(
    title = "Mean |Z| + Variance(|Z|) by Diagnosis",
    x = "Diagnosis",
    y = "L1 Deviation + Variance(|Z|)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

print(p_combined)





# Small epsilon to avoid log(0)
eps <- 1e-12

# Log of mean absolute deviation
scored$log_term1_L1 <- log(scored$term1_L1 + eps)

# Log of variance term
scored$log_term3 <- log(scored$term3_var_absZ + eps)


library(ggplot2)

scored$DIAGNOSIS <- factor(scored$DIAGNOSIS, levels = c("CN","MCI","AD"))
cols <- c("CN"="green", "MCI"="blue", "AD"="red")


p_log_L1 <- ggplot(scored, aes(x = DIAGNOSIS, y = log_term1_L1, fill = DIAGNOSIS)) +
  geom_boxplot(alpha = 0.75, outlier.size = 1.2) +
  scale_fill_manual(values = cols) +
  labs(
    title = "log(Mean Absolute CN-Deviation) by Diagnosis",
    x = "Diagnosis",
    y = "log(Mean |Z|)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

print(p_log_L1)


scored$score_log_sum <- scored$log_term1_L1 + scored$log_term3

p_log_sum <- ggplot(scored, aes(x = DIAGNOSIS, y = score_log_sum, fill = DIAGNOSIS)) +
  geom_boxplot(alpha = 0.75, outlier.size = 1.2) +
  scale_fill_manual(values = cols) +
  labs(
    title = "log(Mean |Z|) + log(Variance(|Z|)) by Diagnosis",
    x = "Diagnosis",
    y = "log(|Z| mean) + log(|Z| variance)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

print(p_log_sum)


## LAST SHOT

# --- Build robust Z-scores vs CN baseline (median/MAD on CN) ---
X <- as.matrix(data_clean[, feature_cols, drop = FALSE])

cn_df <- data_clean[data_clean$DIAGNOSIS == "CN", , drop = FALSE]
m <- apply(cn_df[, feature_cols, drop = FALSE], 2, median, na.rm = TRUE)
s <- apply(cn_df[, feature_cols, drop = FALSE], 2, mad,    na.rm = TRUE)
s[s == 0] <- 1e-6

Z <- sweep(X, 2, m, "-")
Z <- sweep(Z, 2, s, "/")

absZ <- abs(Z)

# Ensure diagnosis ordering/colors for plots
scored$DIAGNOSIS <- factor(scored$DIAGNOSIS, levels = c("CN","MCI","AD"))
cols <- c("CN"="green", "MCI"="blue", "AD"="red")

library(e1071)
library(ggplot2)

# Term 4_1: skewness of |Z| across features for each subject
term4_1_skew_absZ <- apply(absZ, 1, function(v) e1071::skewness(v, na.rm = TRUE, type = 2))
scored$term4_1_skew_absZ <- term4_1_skew_absZ

# Boxplot
ggplot(scored, aes(x = DIAGNOSIS, y = term4_1_skew_absZ, fill = DIAGNOSIS)) +
  geom_boxplot(alpha = 0.75, outlier.size = 1.2) +
  scale_fill_manual(values = cols) +
  labs(
    title = "Term 4_1: Skewness of |Z| by Diagnosis",
    x = "Diagnosis",
    y = "Skewness(|Z|)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")




library(ggplot2)

# Install/load shrinkage covariance helper (recommended)
if (!requireNamespace("corpcor", quietly = TRUE)) install.packages("corpcor")
library(corpcor)

# CN-only Z matrix
Z_cn <- Z[data_clean$DIAGNOSIS == "CN", , drop = FALSE]

# Shrinkage covariance on CN deviations
Sigma_cn <- corpcor::cov.shrink(Z_cn)

# Invert (solve is safe here because cov.shrink is PD)
Sigma_inv <- solve(Sigma_cn)

# Term 4_2: quadratic form z' ??^{-1} z  (per subject)
term4_2_quad <- apply(Z, 1, function(z) as.numeric(t(z) %*% Sigma_inv %*% z))
scored$term4_2_quad <- term4_2_quad

# Boxplot (often needs log scale because it can be heavy-tailed)
ggplot(scored, aes(x = DIAGNOSIS, y = term4_2_quad, fill = DIAGNOSIS)) +
  geom_boxplot(alpha = 0.75, outlier.size = 1.2) +
  scale_fill_manual(values = cols) +
  scale_y_log10() +
  labs(
    title = "Term 4_2: CN Correlation-Pattern Deviation (log scale)",
    x = "Diagnosis",
    y = expression(log[10](z^T*Sigma[CN]^{-1}*z))
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")




## TERM 5

library(stringr)

# Helper: extract region name (ignore ST code suffix)
region_name <- function(x) {
  x %>%
    str_remove("^(Left|Right)\\.") %>%
    str_remove("\\.Volume.*$")
}

# Feature columns only
feature_cols <- setdiff(colnames(data_clean),
                        c("RID","DIAGNOSIS","GENOTYPE","DIAGNOSIS_NUM"))

lr_cols <- feature_cols[grepl("^(Left|Right)\\.", feature_cols)]

lr_df <- data.frame(
  col    = lr_cols,
  side   = ifelse(grepl("^Left\\.", lr_cols), "Left", "Right"),
  region = sapply(lr_cols, region_name),
  stringsAsFactors = FALSE
)

# Build explicit Left/Right pairing table (inner join on region)
left_df  <- lr_df[lr_df$side == "Left",  c("region", "col")]
right_df <- lr_df[lr_df$side == "Right", c("region", "col")]
names(left_df)[2]  <- "left_col"
names(right_df)[2] <- "right_col"

pairs <- merge(left_df, right_df, by = "region")   # keeps only regions that exist on both sides
pairs <- pairs[order(pairs$region), ]              # enforce identical ordering

cat("Number of paired regions:", nrow(pairs), "\n")
cat("Example pairs:\n")
print(head(pairs, 10))

# Extract Z columns in perfectly matched order
Z_left  <- Z[, pairs$left_col,  drop = FALSE]
Z_right <- Z[, pairs$right_col, drop = FALSE]

stopifnot(ncol(Z_left) == ncol(Z_right))

# Term 5: mean absolute hemispheric difference
term5_asym <- rowMeans(abs(Z_left - Z_right), na.rm = TRUE)
scored$term5_asym <- term5_asym




library(ggplot2)

cols <- c("CN" = "green",
          "MCI" = "blue",
          "AD" = "red")

ggplot(scored, aes(x = DIAGNOSIS,
                   y = term5_asym,
                   fill = DIAGNOSIS)) +
  geom_boxplot(alpha = 0.75, outlier.size = 1.2) +
  scale_fill_manual(values = cols) +
  labs(
    title = "Term 5: Hemispheric Asymmetry",
    x = "Diagnosis",
    y = "Mean |Z_left ??? Z_right|"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")

ggplot(scored, aes(x = DIAGNOSIS,
                   y = term5_asym,
                   fill = DIAGNOSIS)) +
  geom_boxplot(alpha = 0.75, outlier.size = 1.2) +
  scale_fill_manual(values = cols) +
  scale_y_log10() +
  labs(
    title = "Term 5: Hemispheric Asymmetry (log scale)",
    x = "Diagnosis",
    y = "Mean |Z_left ??? Z_right| (log scale)"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")



# Make sure all terms exist
terms_df <- scored[, c("log_term1_L1",
                       "log_term3",
                       "term4_1_skew_absZ",
                       "term4_2_quad",
                       "term5_asym")]

cor_mat <- cor(terms_df, use = "complete.obs")

print(round(cor_mat, 4))


cor_t1_t5 <- cor(scored$log_term1_L1,
                 scored$term5_asym,
                 use = "complete.obs")

cat("Correlation(log_term1_L1, term5_asym):",
    round(cor_t1_t5, 4), "\n")

# Select terms
pca_terms <- scored[, c("log_term1_L1",
                        "log_term3",
                        "term4_1_skew_absZ",
                        "term4_2_quad",
                        "term5_asym")]

# Standardize automatically (important!)
pca_fit <- prcomp(pca_terms,
                  center = TRUE,
                  scale. = TRUE)

# Proportion of variance explained
summary(pca_fit)

# Loadings
loadings <- pca_fit$rotation
print(round(loadings[, 1], 4))   # PC1 loadings only

scored$PC1 <- pca_fit$x[, 1]

library(ggplot2)

cols <- c("CN"="green",
          "MCI"="blue",
          "AD"="red")

ggplot(scored, aes(x = DIAGNOSIS,
                   y = PC1,
                   fill = DIAGNOSIS)) +
  geom_boxplot(alpha = 0.75, outlier.size = 1.2) +
  scale_fill_manual(values = cols) +
  labs(
    title = "PC1 of Structural Terms by Diagnosis",
    x = "Diagnosis",
    y = "Principal Component 1"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")


scored$PC2 <- pca_fit$x[, 2]

ggplot(scored, aes(x = PC1, y = PC2, color = DIAGNOSIS)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = cols) +
  theme_bw(base_size = 14) +
  labs(title = "PC1 vs PC2",
       x = "PC1",
       y = "PC2")






library(nnet)     # multinom
library(dplyr)

# Make sure DIAGNOSIS is factor
scored$DIAGNOSIS <- factor(scored$DIAGNOSIS,
                           levels = c("CN","MCI","AD"))

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

cv_multinom_PC1 <- function(df, K = 5, seed = 123) {
  
  folds <- make_stratified_folds(df$DIAGNOSIS, K = K, seed = seed)
  
  acc_vec <- numeric(K)
  all_pred <- rep(NA, nrow(df))
  all_true <- df$DIAGNOSIS
  
  for (fold_id in seq_len(K)) {
    
    test_idx  <- folds[[fold_id]]
    train_idx <- setdiff(seq_len(nrow(df)), test_idx)
    
    trainD <- df[train_idx, ]
    testD  <- df[test_idx, ]
    
    model <- multinom(DIAGNOSIS ~ PC1, data = trainD, trace = FALSE)
    pred  <- predict(model, newdata = testD)
    
    acc_vec[fold_id] <- mean(pred == testD$DIAGNOSIS)
    all_pred[test_idx] <- as.character(pred)
  }
  
  all_pred <- factor(all_pred, levels = levels(all_true))
  
  list(
    mean_accuracy = mean(acc_vec),
    sd_accuracy   = sd(acc_vec),
    confusion     = table(Predicted = all_pred, True = all_true)
  )
}

pc1_res <- cv_multinom_PC1(scored, K = 5)

cat("PC1 multinomial CV mean accuracy:", pc1_res$mean_accuracy, "\n")
cat("PC1 multinomial CV SD:", pc1_res$sd_accuracy, "\n\n")
print(pc1_res$confusion)

