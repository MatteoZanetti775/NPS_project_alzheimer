#---- 0) Load cleaned data ----
data <- read.csv("data.csv")

colnames(data)

## ============================================
## Longitudinal diagnostic plots (2 plots)
## - choose class by code: 1=CN, 2=MCI, 4=AD
## - choose brain region by index within brain columns (start at 17)
## ============================================

library(ggplot2)

## ---- USER CHOICES ----
class_code <- 1   # 1 = CN, 2 = MCI, 4 = AD
r_idx      <- 1   # region index within brain-volume block (1 means column 17)

## ---- SETUP ----
# robust conversion (in case DIAGNOSIS is factor/character)
data$DIAGNOSIS <- as.numeric(as.character(data$DIAGNOSIS))

# brain region columns start at 17
region_start <- 17
region_cols  <- region_start:ncol(data)

# pick region column by index-within-region-block
stopifnot(r_idx >= 1, r_idx <= length(region_cols))
region_col  <- region_cols[r_idx]
region_name <- colnames(data)[region_col]

## ---- FILTER DATA FOR PLOTTING ----
df_plot <- data[
  data$DIAGNOSIS == class_code &
    !is.na(data$month) &
    !is.na(data[[region_col]]),
  c("RID", "month", region_name)   # <-- all by name (no mixing types)
]

## ---- SANITY PRINTS ----
cat("Rows in df_plot:", nrow(df_plot), "\n")
cat("Class code:", class_code, "\n")
cat("Region index within block (r_idx):", r_idx, "\n")
cat("Region column in data (region_col):", region_col, "\n")
cat("Region name:", region_name, "\n")
cat("Diagnosis counts (all data):\n")
print(table(data$DIAGNOSIS, useNA = "ifany"))

if (nrow(df_plot) == 0) {
  stop("df_plot is empty. Check class_code, month, or selected region for missingness.")
}

## ---- PLOT 1: pooled scatter ----
p1 <- ggplot(df_plot, aes(x = month, y = .data[[region_name]])) +
  geom_point(alpha = 0.4) +
  labs(
    title = paste("Class", class_code, "-", region_name),
    x = "Months since baseline",
    y = "Volume"
  ) +
  theme_minimal()

print(p1)

## ---- PLOT 2: spaghetti (per patient) ----
p2 <- ggplot(df_plot, aes(
  x = month,
  y = .data[[region_name]],
  group = RID
)) +
  geom_line(alpha = 0.15) +
  geom_point(alpha = 0.15) +
  labs(
    title = paste("Class", class_code, "-", region_name, "(per patient)"),
    x = "Months since baseline",
    y = "Volume"
  ) +
  theme_minimal()

print(p2)

## ============================================
## Baseline-normalized spaghetti plot
## + mean (red) and variability (dashed: mean ± 1 SD)
## ============================================

library(ggplot2)
library(dplyr)

## ---- USER CHOICES ----
class_code <- 4   # 1 = CN, 2 = MCI, 4 = AD
r_idx      <- 1   # region index within brain columns

## ---- SETUP ----
data$DIAGNOSIS <- as.numeric(as.character(data$DIAGNOSIS))

region_start <- 17
region_cols  <- region_start:ncol(data)
region_col   <- region_cols[r_idx]
region_name  <- colnames(data)[region_col]

## ---- FILTER ----
df_plot <- data[
  data$DIAGNOSIS == class_code &
    !is.na(data$month) &
    !is.na(data[[region_col]]),
  c("RID", "month", region_name)
]

## ---- BASELINE NORMALIZATION (percent change) ----
df_plot <- df_plot %>%
  arrange(RID, month) %>%
  group_by(RID) %>%
  mutate(
    baseline   = first(.data[[region_name]]),
    pct_change = 100 * (.data[[region_name]] - baseline) / baseline
  ) %>%
  ungroup()

## ---- MEAN + VARIABILITY AT EACH MONTH ----
df_summary <- df_plot %>%
  group_by(month) %>%
  summarise(
    n = sum(!is.na(pct_change)),
    mean_pct = mean(pct_change, na.rm = TRUE),
    sd_pct   = sd(pct_change, na.rm = TRUE),
    # variance if you want it later:
    var_pct  = var(pct_change, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(month) %>%
  mutate(
    upper_1sd = mean_pct + sd_pct,
    lower_1sd = mean_pct - sd_pct
  )

## ---- PLOT: spaghetti + mean + dashed ±1SD ----
ggplot(df_plot, aes(x = month, y = pct_change, group = RID)) +
  geom_line(alpha = 0.15) +
  geom_point(alpha = 0.15) +
  
  # Mean (red, thicker)
  geom_line(
    data = df_summary,
    aes(x = month, y = mean_pct),
    color = "red",
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  
  # Variability bands (dashed): mean ± 1 SD
  geom_line(
    data = df_summary,
    aes(x = month, y = upper_1sd),
    linetype = "dashed",
    linewidth = 1.0,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = df_summary,
    aes(x = month, y = lower_1sd),
    linetype = "dashed",
    linewidth = 1.0,
    inherit.aes = FALSE
  ) +
  
  labs(
    title = paste("Class", class_code, "-", region_name, "(% change: spaghetti + mean ± 1 SD)"),
    x = "Months since baseline",
    y = "Percent change"
  ) +
  theme_minimal()

## ============================================
## Spaghetti MEAN comparison: CN vs MCI vs AD
## (baseline-normalized, no smoothing)
## ============================================

library(ggplot2)
library(dplyr)

## ---- USER CHOICES ----
r_idx <- 4   # region index within brain columns (1 = column 17)

## ---- SETUP ----
data$DIAGNOSIS <- as.numeric(as.character(data$DIAGNOSIS))

region_start <- 17
region_cols  <- region_start:ncol(data)
region_col   <- region_cols[r_idx]
region_name  <- colnames(data)[region_col]

## ---- KEEP ONLY CLASSES OF INTEREST ----
class_codes <- c(1, 2, 4)  # CN, MCI, AD

df_plot <- data[
  data$DIAGNOSIS %in% class_codes &
    !is.na(data$month) &
    !is.na(data[[region_col]]),
  c("RID", "month", "DIAGNOSIS", region_name)
]

## ---- BASELINE NORMALIZATION (percent change) ----
df_plot <- df_plot %>%
  arrange(RID, month) %>%
  group_by(RID) %>%
  mutate(
    baseline   = first(.data[[region_name]]),
    pct_change = 100 * (.data[[region_name]] - baseline) / baseline
  ) %>%
  ungroup()

## ---- COMPUTE SPAGHETTI MEAN PER CLASS ----
df_mean <- df_plot %>%
  group_by(DIAGNOSIS, month) %>%
  summarise(
    mean_pct = mean(pct_change, na.rm = TRUE),
    n = sum(!is.na(pct_change)),
    .groups = "drop"
  ) %>%
  arrange(DIAGNOSIS, month)

## ---- LABELS FOR PLOTTING ----
df_mean$CLASS <- factor(
  df_mean$DIAGNOSIS,
  levels = c(1, 2, 4),
  labels = c("CN", "MCI", "AD")
)

## ---- PLOT: three spaghetti means together ----
ggplot(df_mean, aes(
  x = month,
  y = mean_pct,
  color = CLASS
)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2) +
  labs(
    title = paste("Spaghetti mean comparison (% change) -", region_name),
    x = "Months since baseline",
    y = "Mean percent change from baseline",
    color = "Class"
  ) +
  theme_minimal()






library(dplyr)
library(mgcv)

set.seed(123)

K_folds     <- 5
max_visits  <- 4
eps_var     <- 1e-4
min_points  <- 30
k_spline    <- 10

region_start <- 17
region_cols  <- region_start:ncol(data)
 region_cols <- region_cols[1:30]  # optional speed-up

classes <- c(1, 2, 4)

data$DIAGNOSIS <- as.numeric(as.character(data$DIAGNOSIS))
data$month     <- as.numeric(as.character(data$month))

keep_cols <- c("RID", "month", "DIAGNOSIS", colnames(data)[region_cols])
df0 <- data[, keep_cols]

df0 <- df0[!is.na(df0$RID) & !is.na(df0$month) & !is.na(df0$DIAGNOSIS), ]
df0 <- df0[df0$DIAGNOSIS %in% classes, ]

# IMPORTANT: ensure column names are unique (prevents silent weirdness)
colnames(df0) <- make.unique(colnames(df0))

# define region names ONCE from df0 (the object we actually use)
region_names <- colnames(df0)[4:ncol(df0)]  # after RID,month,DIAGNOSIS

label_class <- function(x) factor(x, levels = c(1,2,4), labels = c("CN","MCI","AD"))

all_rids <- sort(unique(df0$RID))
fold_id  <- sample(rep(1:K_folds, length.out = length(all_rids)))
rid2fold <- data.frame(RID = all_rids, fold = fold_id)

all_preds <- list()

for (fold in 1:K_folds) {
  cat("\n=============================\n")
  cat("Fold", fold, "of", K_folds, "\n")
  cat("=============================\n")
  
  test_rids  <- rid2fold$RID[rid2fold$fold == fold]
  train_rids <- rid2fold$RID[rid2fold$fold != fold]
  
  train_df <- df0[df0$RID %in% train_rids, ]
  test_df  <- df0[df0$RID %in% test_rids, ]
  
  # models[[region_name]][[class]] = list(mean=..., logvar=...)
  models <- setNames(vector("list", length(region_names)), region_names)
  
  pb_reg <- txtProgressBar(min = 0, max = length(region_names), style = 3)
  
  for (ri in seq_along(region_names)) {
    region_name <- region_names[ri]
    
    # extra guard (shouldn't trigger now, but safe)
    if (!region_name %in% names(train_df)) {
      models[[region_name]] <- NULL
      setTxtProgressBar(pb_reg, ri)
      next
    }
    
    tmp <- train_df[, c("RID","month","DIAGNOSIS", region_name)]
    colnames(tmp)[4] <- "y"
    
    tmp <- tmp[!is.na(tmp$y), ]
    if (nrow(tmp) < min_points) {
      models[[region_name]] <- NULL
      setTxtProgressBar(pb_reg, ri)
      next
    }
    
    tmp <- tmp %>%
      arrange(RID, month) %>%
      group_by(RID) %>%
      mutate(
        baseline = first(y),
        pct_change = 100 * (y - baseline) / baseline
      ) %>%
      ungroup()
    
    reg_models <- list()
    
    for (k in classes) {
      subk <- tmp[tmp$DIAGNOSIS == k, ]
      if (nrow(subk) < min_points) {
        reg_models[[as.character(k)]] <- NULL
        next
      }
      
      m_mean <- tryCatch(
        gam(pct_change ~ s(month, k = k_spline), data = subk, method = "REML"),
        error = function(e) NULL
      )
      if (is.null(m_mean)) {
        reg_models[[as.character(k)]] <- NULL
        next
      }
      
      resid <- subk$pct_change - predict(m_mean, newdata = subk)
      subk$log_r2 <- log(pmax(resid^2, 1e-8))
      
      m_logvar <- tryCatch(
        gam(log_r2 ~ s(month, k = k_spline), data = subk, method = "REML"),
        error = function(e) NULL
      )
      if (is.null(m_logvar)) {
        reg_models[[as.character(k)]] <- NULL
        next
      }
      
      reg_models[[as.character(k)]] <- list(mean = m_mean, logvar = m_logvar)
    }
    
    models[[region_name]] <- reg_models
    setTxtProgressBar(pb_reg, ri)
  }
  close(pb_reg)
  
  # ---- prediction loop ----
  test_patients <- sort(unique(test_df$RID))
  pb_pat <- txtProgressBar(min = 0, max = length(test_patients), style = 3)
  
  fold_out <- data.frame(RID = integer(0), true = integer(0), pred = integer(0))
  
  train_majority <- as.numeric(names(sort(table(train_df$DIAGNOSIS), decreasing = TRUE))[1])
  
  for (pi in seq_along(test_patients)) {
    rid <- test_patients[pi]
    subp <- test_df[test_df$RID == rid, ]
    subp <- subp[order(subp$month), ]
    
    true_class <- subp$DIAGNOSIS[1]
    subp <- subp[1:min(nrow(subp), max_visits), ]
    
    D <- setNames(rep(0, length(classes)), as.character(classes))
    used_any <- FALSE
    
    for (region_name in region_names) {
      reg_models <- models[[region_name]]
      if (is.null(reg_models)) next
      
      yvec <- subp[[region_name]]
      if (all(is.na(yvec))) next
      
      idx_ok <- which(!is.na(yvec))
      if (length(idx_ok) < 2) next
      
      tvec <- subp$month[idx_ok]
      yok  <- yvec[idx_ok]
      base <- yok[1]
      zvec <- 100 * (yok - base) / base
      
      for (k in classes) {
        mk <- reg_models[[as.character(k)]]
        if (is.null(mk)) next
        
        mu_hat <- predict(mk$mean,   newdata = data.frame(month = tvec))
        lv_hat <- predict(mk$logvar, newdata = data.frame(month = tvec))
        sig2   <- pmax(exp(lv_hat), eps_var)
        
        D[as.character(k)] <- D[as.character(k)] + sum((zvec - mu_hat)^2 / (sig2 + eps_var))
        used_any <- TRUE
      }
    }
    
    pred_class <- if (!used_any || any(!is.finite(D))) train_majority else as.numeric(names(which.min(D)))
    
    fold_out <- rbind(fold_out, data.frame(RID = rid, true = true_class, pred = pred_class))
    setTxtProgressBar(pb_pat, pi)
  }
  close(pb_pat)
  
  all_preds[[fold]] <- fold_out
  
  cat("\nFold accuracy:", round(mean(fold_out$true == fold_out$pred), 4), "\n")
  cat("Confusion (rows=true, cols=pred):\n")
  print(table(True = label_class(fold_out$true), Pred = label_class(fold_out$pred)))
}

preds <- do.call(rbind, all_preds)

cat("\n=============================\n")
cat("Overall CV accuracy:", round(mean(preds$true == preds$pred), 4), "\n")
cat("=============================\n")
cat("Overall confusion (rows=true, cols=pred):\n")
print(table(True = label_class(preds$true), Pred = label_class(preds$pred)))









## ============================================================
## PURELY FUNCTIONAL 3-class classification via Functional Depth
## Fraiman-Muniz depth (integrated pointwise rank depth)
##
## - 3 classes: 1=CN, 2=MCI, 4=AD
## - Uses <=4 earliest visits per patient (timely prediction)
## - Represents each patient x region as a function via interpolation
## - Classify new patient by max depth-to-class (summed over regions)
## - Patient-level K-fold CV + progress bars
##
## Requirements: dplyr
## ============================================================

library(dplyr)

set.seed(123)

## -------------------------
## USER SETTINGS
## -------------------------
K_folds      <- 5
max_visits   <- 4
classes      <- c(1, 2, 4)

region_start <- 17
region_cols  <- region_start:ncol(data)

# Optional speed-up for first run:
region_cols <- region_cols[1:30]

# Common time grid (months). Choose a step that's not too fine.
grid_step    <- 6

# Minimum number of training curves required at a grid timepoint to score depth
min_m_at_t   <- 10

# Whether to baseline-normalize (purely functional operator on curves).
# TRUE sets each patient's curve to percent-change-from-first-observed-value.
use_pct_change <- TRUE

## -------------------------
## BASIC CLEANUP
## -------------------------
data$DIAGNOSIS <- as.numeric(as.character(data$DIAGNOSIS))
data$month     <- as.numeric(as.character(data$month))

# Keep only needed columns
keep_cols <- c("RID", "month", "DIAGNOSIS", colnames(data)[region_cols])
df0 <- data[, keep_cols]
df0 <- df0[!is.na(df0$RID) & !is.na(df0$month) & !is.na(df0$DIAGNOSIS), ]
df0 <- df0[df0$DIAGNOSIS %in% classes, ]

# Ensure unique names (robust)
colnames(df0) <- make.unique(colnames(df0))

region_names <- colnames(df0)[4:ncol(df0)]

label_class <- function(x) factor(x, levels = c(1,2,4), labels = c("CN","MCI","AD"))

## -------------------------
## Create patient folds
## -------------------------
all_rids <- sort(unique(df0$RID))
fold_id  <- sample(rep(1:K_folds, length.out = length(all_rids)))
rid2fold <- data.frame(RID = all_rids, fold = fold_id)

## ============================================================
## Helper: build <=4-visit trajectories per patient for one region
## Returns a list: curves_by_class[[k]] is a matrix (n_pat x n_grid)
## ============================================================

# linear interpolation onto grid; outside range -> NA
interp_to_grid <- function(t, y, grid) {
  if (length(t) < 2) return(rep(NA_real_, length(grid)))
  approx(x = t, y = y, xout = grid, method = "linear", rule = 1, ties = "ordered")$y
}

# Fraiman-Muniz pointwise depth at a single time t:
# depth = 1 - 2*|F(y)-0.5| where F is empirical CDF from training values
fm_point_depth <- function(y0, y_train_vec) {
  y_train_vec <- y_train_vec[is.finite(y_train_vec)]
  m <- length(y_train_vec)
  if (!is.finite(y0) || m < 2) return(NA_real_)
  # empirical CDF value at y0 (using <=)
  F <- mean(y_train_vec <= y0)
  1 - 2 * abs(F - 0.5)
}

# FM functional depth of one curve vs training curves (matrix: n_train x n_grid)
fm_depth_curve <- function(y0_grid, Y_train_mat, min_m_at_t = 10) {
  n_grid <- length(y0_grid)
  d_t <- rep(NA_real_, n_grid)
  
  for (g in seq_len(n_grid)) {
    train_col <- Y_train_mat[, g]
    train_col <- train_col[is.finite(train_col)]
    if (length(train_col) < min_m_at_t) next
    d_t[g] <- fm_point_depth(y0_grid[g], train_col)
  }
  
  if (all(is.na(d_t))) return(NA_real_)
  mean(d_t, na.rm = TRUE)
}

## ============================================================
## CROSS-VALIDATION
## ============================================================
all_preds <- list()

for (fold in 1:K_folds) {
  cat("\n=============================\n")
  cat("Fold", fold, "of", K_folds, "\n")
  cat("=============================\n")
  
  test_rids  <- rid2fold$RID[rid2fold$fold == fold]
  train_rids <- rid2fold$RID[rid2fold$fold != fold]
  
  train_df <- df0[df0$RID %in% train_rids, ]
  test_df  <- df0[df0$RID %in% test_rids, ]
  
  ## ---- Determine a sensible grid based on <=4 visits only (train set) ----
  # Use each patient's first <=4 visits to compute max month to cover
  train_first4 <- train_df %>%
    arrange(RID, month) %>%
    group_by(RID) %>%
    slice_head(n = max_visits) %>%
    ungroup()
  
  Tmax <- max(train_first4$month, na.rm = TRUE)
  grid <- seq(0, Tmax, by = grid_step)
  
  cat("Grid: 0 to", Tmax, "by", grid_step, "months (n =", length(grid), ")\n")
  
  ## ============================================================
  ## Build training curve matrices per region and class
  ## curves_train[[region]][[class]] = matrix(n_pat_in_class x n_grid)
  ## Also store patient order for each class for debug (optional)
  ## ============================================================
  curves_train <- setNames(vector("list", length(region_names)), region_names)
  
  pb_reg <- txtProgressBar(min = 0, max = length(region_names), style = 3)
  
  for (ri in seq_along(region_names)) {
    region <- region_names[ri]
    
    # Use <=4 visits only for training curves
    tmp <- train_df[, c("RID","month","DIAGNOSIS", region)]
    colnames(tmp)[4] <- "y"
    tmp <- tmp[!is.na(tmp$y), ]
    
    # Keep <=4 earliest per RID
    tmp <- tmp %>%
      arrange(RID, month) %>%
      group_by(RID) %>%
      slice_head(n = max_visits) %>%
      ungroup()
    
    if (nrow(tmp) == 0) {
      curves_train[[region]] <- NULL
      setTxtProgressBar(pb_reg, ri)
      next
    }
    
    # Build curve matrices per class
    reg_by_class <- list()
    
    for (k in classes) {
      subk <- tmp[tmp$DIAGNOSIS == k, ]
      rids_k <- sort(unique(subk$RID))
      if (length(rids_k) < 5) {
        reg_by_class[[as.character(k)]] <- NULL
        next
      }
      
      Yk <- matrix(NA_real_, nrow = length(rids_k), ncol = length(grid))
      
      for (ii in seq_along(rids_k)) {
        rid <- rids_k[ii]
        si <- subk[subk$RID == rid, ]
        si <- si[order(si$month), ]
        t <- si$month
        y <- si$y
        
        # baseline-normalize if requested (percent change from first observed value)
        if (use_pct_change) {
          if (length(y) >= 1 && is.finite(y[1]) && y[1] != 0) {
            y <- 100 * (y - y[1]) / y[1]
          } else {
            # if baseline invalid, skip
            next
          }
        }
        
        Yk[ii, ] <- interp_to_grid(t, y, grid)
      }
      
      reg_by_class[[as.character(k)]] <- Yk
    }
    
    curves_train[[region]] <- reg_by_class
    setTxtProgressBar(pb_reg, ri)
  }
  close(pb_reg)
  
  ## ============================================================
  ## Predict test patients
  ## ============================================================
  test_patients <- sort(unique(test_df$RID))
  pb_pat <- txtProgressBar(min = 0, max = length(test_patients), style = 3)
  
  fold_out <- data.frame(RID = integer(0), true = integer(0), pred = integer(0))
  
  # fallback: majority class in training (safety)
  majority_train <- as.numeric(names(sort(table(train_df$DIAGNOSIS), decreasing = TRUE))[1])
  
  for (pi in seq_along(test_patients)) {
    rid <- test_patients[pi]
    subp <- test_df[test_df$RID == rid, ]
    subp <- subp[order(subp$month), ]
    
    true_class <- subp$DIAGNOSIS[1]
    subp <- subp[1:min(nrow(subp), max_visits), ]
    
    # total depth score per class (sum across regions)
    S <- setNames(rep(0, length(classes)), as.character(classes))
    used_any <- FALSE
    
    for (region in region_names) {
      reg_models <- curves_train[[region]]
      if (is.null(reg_models)) next
      
      # patient y in this region
      yvec <- subp[[region]]
      if (all(is.na(yvec))) next
      
      idx_ok <- which(!is.na(yvec))
      if (length(idx_ok) < 2) next
      
      t <- subp$month[idx_ok]
      y <- yvec[idx_ok]
      
      if (use_pct_change) {
        if (!is.finite(y[1]) || y[1] == 0) next
        y <- 100 * (y - y[1]) / y[1]
      }
      
      y_grid <- interp_to_grid(t, y, grid)
      
      # compute depth against each class for THIS region
      for (k in classes) {
        Y_train <- reg_models[[as.character(k)]]
        if (is.null(Y_train)) next
        
        d_fk <- fm_depth_curve(y_grid, Y_train, min_m_at_t = min_m_at_t)
        if (is.finite(d_fk)) {
          S[as.character(k)] <- S[as.character(k)] + d_fk
          used_any <- TRUE
        }
      }
    }
    
    pred_class <- if (!used_any || any(!is.finite(S))) {
      majority_train
    } else {
      as.numeric(names(which.max(S)))
    }
    
    fold_out <- rbind(fold_out, data.frame(RID = rid, true = true_class, pred = pred_class))
    setTxtProgressBar(pb_pat, pi)
  }
  close(pb_pat)
  
  all_preds[[fold]] <- fold_out
  
  cat("\nFold accuracy:", round(mean(fold_out$true == fold_out$pred), 4), "\n")
  cat("Confusion (rows=true, cols=pred):\n")
  print(table(True = label_class(fold_out$true), Pred = label_class(fold_out$pred)))
}

## ============================================================
## Overall performance
## ============================================================
preds <- do.call(rbind, all_preds)

cat("\n=============================\n")
cat("Overall CV accuracy:", round(mean(preds$true == preds$pred), 4), "\n")
cat("=============================\n")
cat("Overall confusion (rows=true, cols=pred):\n")
print(table(True = label_class(preds$true), Pred = label_class(preds$pred)))











library(dplyr)
library(class)

set.seed(123)

K_folds        <- 5
max_visits     <- 4
classes        <- c(1,2,4)

region_start   <- 17
region_cols    <- region_start:ncol(data)
# region_cols <- region_cols[1:30]  # optional speed-up

grid_step      <- 6
Tmax_use       <- 48     # <<< IMPORTANT: early horizon cap
min_m_at_t     <- 10
use_pct_change <- TRUE

k_nn           <- 15     # k in kNN in depth-space

data$DIAGNOSIS <- as.numeric(as.character(data$DIAGNOSIS))
data$month     <- as.numeric(as.character(data$month))

keep_cols <- c("RID","month","DIAGNOSIS", colnames(data)[region_cols])
df0 <- data[, keep_cols]
df0 <- df0[!is.na(df0$RID) & !is.na(df0$month) & !is.na(df0$DIAGNOSIS), ]
df0 <- df0[df0$DIAGNOSIS %in% classes, ]
colnames(df0) <- make.unique(colnames(df0))
region_names <- colnames(df0)[4:ncol(df0)]

label_class <- function(x) factor(x, levels=c(1,2,4), labels=c("CN","MCI","AD"))

all_rids <- sort(unique(df0$RID))
fold_id  <- sample(rep(1:K_folds, length.out = length(all_rids)))
rid2fold <- data.frame(RID = all_rids, fold = fold_id)

interp_to_grid <- function(t, y, grid) {
  if (length(t) < 2) return(rep(NA_real_, length(grid)))
  approx(x=t, y=y, xout=grid, method="linear", rule=1, ties="ordered")$y
}

fm_point_depth <- function(y0, y_train_vec) {
  y_train_vec <- y_train_vec[is.finite(y_train_vec)]
  m <- length(y_train_vec)
  if (!is.finite(y0) || m < 2) return(NA_real_)
  F <- mean(y_train_vec <= y0)
  1 - 2 * abs(F - 0.5)
}

fm_depth_curve <- function(y0_grid, Y_train_mat, min_m_at_t = 10) {
  n_grid <- length(y0_grid)
  d_t <- rep(NA_real_, n_grid)
  for (g in seq_len(n_grid)) {
    train_col <- Y_train_mat[, g]
    train_col <- train_col[is.finite(train_col)]
    if (length(train_col) < min_m_at_t) next
    d_t[g] <- fm_point_depth(y0_grid[g], train_col)
  }
  if (all(is.na(d_t))) return(NA_real_)
  mean(d_t, na.rm=TRUE)
}

all_preds <- list()

for (fold in 1:K_folds) {
  cat("\n=============================\n")
  cat("Fold", fold, "of", K_folds, "\n")
  cat("=============================\n")
  
  test_rids  <- rid2fold$RID[rid2fold$fold == fold]
  train_rids <- rid2fold$RID[rid2fold$fold != fold]
  
  train_df <- df0[df0$RID %in% train_rids, ]
  test_df  <- df0[df0$RID %in% test_rids, ]
  
  grid <- seq(0, Tmax_use, by = grid_step)
  
  ## --- Build training curves per region/class using <=4 visits ---
  curves_train <- setNames(vector("list", length(region_names)), region_names)
  
  pb_reg <- txtProgressBar(min=0, max=length(region_names), style=3)
  for (ri in seq_along(region_names)) {
    region <- region_names[ri]
    tmp <- train_df[, c("RID","month","DIAGNOSIS", region)]
    colnames(tmp)[4] <- "y"
    tmp <- tmp[!is.na(tmp$y), ]
    
    tmp <- tmp %>%
      arrange(RID, month) %>%
      group_by(RID) %>%
      slice_head(n = max_visits) %>%
      ungroup()
    
    reg_by_class <- list()
    for (k in classes) {
      subk <- tmp[tmp$DIAGNOSIS == k, ]
      rids_k <- sort(unique(subk$RID))
      if (length(rids_k) < 5) { reg_by_class[[as.character(k)]] <- NULL; next }
      
      Yk <- matrix(NA_real_, nrow=length(rids_k), ncol=length(grid))
      for (ii in seq_along(rids_k)) {
        rid <- rids_k[ii]
        si <- subk[subk$RID == rid, ]
        si <- si[order(si$month), ]
        t <- si$month; y <- si$y
        
        if (use_pct_change) {
          if (!is.finite(y[1]) || y[1] == 0) next
          y <- 100 * (y - y[1]) / y[1]
        }
        
        Yk[ii, ] <- interp_to_grid(t, y, grid)
      }
      reg_by_class[[as.character(k)]] <- Yk
    }
    
    curves_train[[region]] <- reg_by_class
    setTxtProgressBar(pb_reg, ri)
  }
  close(pb_reg)
  
  ## --- Build depth-vectors for TRAIN patients (DD features) ---
  train_patients <- sort(unique(train_df$RID))
  X_train <- matrix(NA_real_, nrow=length(train_patients), ncol=length(classes))
  colnames(X_train) <- c("d_CN","d_MCI","d_AD")
  y_train <- integer(length(train_patients))
  
  pb_tr <- txtProgressBar(min=0, max=length(train_patients), style=3)
  for (i in seq_along(train_patients)) {
    rid <- train_patients[i]
    subp <- train_df[train_df$RID == rid, ]
    subp <- subp[order(subp$month), ]
    y_train[i] <- subp$DIAGNOSIS[1]
    subp <- subp[1:min(nrow(subp), max_visits), ]
    
    S <- setNames(rep(0, length(classes)), as.character(classes))
    used <- FALSE
    
    for (region in region_names) {
      reg_models <- curves_train[[region]]
      if (is.null(reg_models)) next
      
      yvec <- subp[[region]]
      if (all(is.na(yvec))) next
      idx_ok <- which(!is.na(yvec))
      if (length(idx_ok) < 2) next
      
      t <- subp$month[idx_ok]
      y <- yvec[idx_ok]
      
      if (use_pct_change) {
        if (!is.finite(y[1]) || y[1] == 0) next
        y <- 100 * (y - y[1]) / y[1]
      }
      
      y_grid <- interp_to_grid(t, y, grid)
      
      for (k in classes) {
        Y_train_k <- reg_models[[as.character(k)]]
        if (is.null(Y_train_k)) next
        d_fk <- fm_depth_curve(y_grid, Y_train_k, min_m_at_t = min_m_at_t)
        if (is.finite(d_fk)) {
          S[as.character(k)] <- S[as.character(k)] + d_fk
          used <- TRUE
        }
      }
    }
    
    if (used) {
      X_train[i, ] <- c(S["1"], S["2"], S["4"])
    }
    setTxtProgressBar(pb_tr, i)
  }
  close(pb_tr)
  
  keep <- complete.cases(X_train) & is.finite(y_train)
  X_train <- X_train[keep, , drop=FALSE]
  y_train <- y_train[keep]
  
  ## --- Build depth-vectors for TEST patients and classify by kNN ---
  test_patients <- sort(unique(test_df$RID))
  pb_te <- txtProgressBar(min=0, max=length(test_patients), style=3)
  
  fold_out <- data.frame(RID=integer(0), true=integer(0), pred=integer(0))
  
  for (i in seq_along(test_patients)) {
    rid <- test_patients[i]
    subp <- test_df[test_df$RID == rid, ]
    subp <- subp[order(subp$month), ]
    true <- subp$DIAGNOSIS[1]
    subp <- subp[1:min(nrow(subp), max_visits), ]
    
    S <- setNames(rep(0, length(classes)), as.character(classes))
    used <- FALSE
    
    for (region in region_names) {
      reg_models <- curves_train[[region]]
      if (is.null(reg_models)) next
      
      yvec <- subp[[region]]
      if (all(is.na(yvec))) next
      idx_ok <- which(!is.na(yvec))
      if (length(idx_ok) < 2) next
      
      t <- subp$month[idx_ok]
      y <- yvec[idx_ok]
      
      if (use_pct_change) {
        if (!is.finite(y[1]) || y[1] == 0) next
        y <- 100 * (y - y[1]) / y[1]
      }
      
      y_grid <- interp_to_grid(t, y, grid)
      
      for (k in classes) {
        Y_train_k <- reg_models[[as.character(k)]]
        if (is.null(Y_train_k)) next
        d_fk <- fm_depth_curve(y_grid, Y_train_k, min_m_at_t = min_m_at_t)
        if (is.finite(d_fk)) {
          S[as.character(k)] <- S[as.character(k)] + d_fk
          used <- TRUE
        }
      }
    }
    
    if (!used) {
      pred <- as.numeric(names(sort(table(train_df$DIAGNOSIS), decreasing=TRUE))[1])
    } else {
      x_test <- matrix(c(S["1"], S["2"], S["4"]), nrow=1)
      colnames(x_test) <- colnames(X_train)
      pred <- knn(train = X_train, test = x_test, cl = y_train, k = k_nn)
      pred <- as.numeric(as.character(pred))
    }
    
    fold_out <- rbind(fold_out, data.frame(RID=rid, true=true, pred=pred))
    setTxtProgressBar(pb_te, i)
  }
  close(pb_te)
  
  all_preds[[fold]] <- fold_out
  cat("\nFold accuracy:", round(mean(fold_out$true == fold_out$pred), 4), "\n")
  cat("Confusion (rows=true, cols=pred):\n")
  print(table(True = label_class(fold_out$true), Pred = label_class(fold_out$pred)))
}

preds <- do.call(rbind, all_preds)
cat("\n=============================\n")
cat("Overall CV accuracy:", round(mean(preds$true == preds$pred), 4), "\n")
cat("=============================\n")
cat("Overall confusion (rows=true, cols=pred):\n")
print(table(True = label_class(preds$true), Pred = label_class(preds$pred)))









## ============================================================
## PURELY FUNCTIONAL 3-class classification via Depth-Depth (DD)
## + Class-conditional depth normalization (Idea #1)
##
## - 3 classes: 1=CN, 2=MCI, 4=AD
## - Uses <=4 earliest visits per patient (timely prediction)
## - Curves built by linear interpolation on a common grid
## - Depth: Fraiman-Muniz (integrated pointwise rank depth)
## - DD-features: (Depth_to_CN, Depth_to_MCI, Depth_to_AD)
## - NORMALIZATION: each depth divided by typical in-class self-depth
##   so depths are comparable across classes
## - Classifier: kNN in depth space
## - Patient-level K-fold CV + progress bars
##
## Packages: dplyr, class
## ============================================================

library(dplyr)
library(class)

set.seed(123)

## -------------------------
## USER SETTINGS
## -------------------------
K_folds        <- 5
max_visits     <- 4
classes        <- c(1, 2, 4)

region_start   <- 17
region_cols    <- region_start:ncol(data)

# Optional speed-up for first run:
# region_cols <- region_cols[1:30]

grid_step      <- 6
Tmax_use       <- 48     # early horizon cap (recommended)
min_m_at_t     <- 10     # min training curves available at timepoint
use_pct_change <- TRUE   # baseline-normalize per patient (percent change)
k_nn           <- 15     # k in kNN in depth-space

## Typical depth estimator: "mean" or "median" (median is more robust)
typical_stat   <- "median"

## Stabilizer to avoid division by zero
eps_typical    <- 1e-6

## -------------------------
## BASIC CLEANUP
## -------------------------
data$DIAGNOSIS <- as.numeric(as.character(data$DIAGNOSIS))
data$month     <- as.numeric(as.character(data$month))

keep_cols <- c("RID","month","DIAGNOSIS", colnames(data)[region_cols])
df0 <- data[, keep_cols]
df0 <- df0[!is.na(df0$RID) & !is.na(df0$month) & !is.na(df0$DIAGNOSIS), ]
df0 <- df0[df0$DIAGNOSIS %in% classes, ]

colnames(df0) <- make.unique(colnames(df0))
region_names <- colnames(df0)[4:ncol(df0)]

label_class <- function(x) factor(x, levels=c(1,2,4), labels=c("CN","MCI","AD"))

## -------------------------
## Create patient folds
## -------------------------
all_rids <- sort(unique(df0$RID))
fold_id  <- sample(rep(1:K_folds, length.out = length(all_rids)))
rid2fold <- data.frame(RID = all_rids, fold = fold_id)

## -------------------------
## Helper functions
## -------------------------
interp_to_grid <- function(t, y, grid) {
  if (length(t) < 2) return(rep(NA_real_, length(grid)))
  approx(x=t, y=y, xout=grid, method="linear", rule=1, ties="ordered")$y
}

fm_point_depth <- function(y0, y_train_vec) {
  y_train_vec <- y_train_vec[is.finite(y_train_vec)]
  m <- length(y_train_vec)
  if (!is.finite(y0) || m < 2) return(NA_real_)
  F <- mean(y_train_vec <= y0)
  1 - 2 * abs(F - 0.5)
}

fm_depth_curve <- function(y0_grid, Y_train_mat, min_m_at_t = 10) {
  n_grid <- length(y0_grid)
  d_t <- rep(NA_real_, n_grid)
  
  for (g in seq_len(n_grid)) {
    train_col <- Y_train_mat[, g]
    train_col <- train_col[is.finite(train_col)]
    if (length(train_col) < min_m_at_t) next
    d_t[g] <- fm_point_depth(y0_grid[g], train_col)
  }
  
  if (all(is.na(d_t))) return(NA_real_)
  mean(d_t, na.rm=TRUE)
}

get_typical <- function(x, stat = "median") {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  if (stat == "mean") return(mean(x))
  median(x)
}

## ============================================================
## CROSS-VALIDATION
## ============================================================
all_preds <- list()

for (fold in 1:K_folds) {
  cat("\n=============================\n")
  cat("Fold", fold, "of", K_folds, "\n")
  cat("=============================\n")
  
  test_rids  <- rid2fold$RID[rid2fold$fold == fold]
  train_rids <- rid2fold$RID[rid2fold$fold != fold]
  
  train_df <- df0[df0$RID %in% train_rids, ]
  test_df  <- df0[df0$RID %in% test_rids, ]
  
  grid <- seq(0, Tmax_use, by = grid_step)
  
  ## ----------------------------------------------------------
  ## Build training curves per region/class using <=4 visits
  ## curves_train[[region]][[class]] = matrix(n_pat_in_class x n_grid)
  ## ----------------------------------------------------------
  curves_train <- setNames(vector("list", length(region_names)), region_names)
  
  pb_reg <- txtProgressBar(min=0, max=length(region_names), style=3)
  for (ri in seq_along(region_names)) {
    region <- region_names[ri]
    
    tmp <- train_df[, c("RID","month","DIAGNOSIS", region)]
    colnames(tmp)[4] <- "y"
    tmp <- tmp[!is.na(tmp$y), ]
    
    tmp <- tmp %>%
      arrange(RID, month) %>%
      group_by(RID) %>%
      slice_head(n = max_visits) %>%
      ungroup()
    
    reg_by_class <- list()
    
    for (k in classes) {
      subk <- tmp[tmp$DIAGNOSIS == k, ]
      rids_k <- sort(unique(subk$RID))
      if (length(rids_k) < 5) { reg_by_class[[as.character(k)]] <- NULL; next }
      
      Yk <- matrix(NA_real_, nrow=length(rids_k), ncol=length(grid))
      
      for (ii in seq_along(rids_k)) {
        rid <- rids_k[ii]
        si <- subk[subk$RID == rid, ]
        si <- si[order(si$month), ]
        t <- si$month
        y <- si$y
        
        if (use_pct_change) {
          if (!is.finite(y[1]) || y[1] == 0) next
          y <- 100 * (y - y[1]) / y[1]
        }
        
        Yk[ii, ] <- interp_to_grid(t, y, grid)
      }
      
      reg_by_class[[as.character(k)]] <- Yk
    }
    
    curves_train[[region]] <- reg_by_class
    setTxtProgressBar(pb_reg, ri)
  }
  close(pb_reg)
  
  ## ----------------------------------------------------------
  ## Compute class-conditional typical self-depths per region
  ## typical_depth[[region]][[class]] = typical depth of class-k curves w.r.t class-k
  ## ----------------------------------------------------------
  typical_depth <- setNames(vector("list", length(region_names)), region_names)
  
  pb_typ <- txtProgressBar(min=0, max=length(region_names), style=3)
  for (ri in seq_along(region_names)) {
    region <- region_names[ri]
    reg_by_class <- curves_train[[region]]
    if (is.null(reg_by_class)) {
      typical_depth[[region]] <- NULL
      setTxtProgressBar(pb_typ, ri)
      next
    }
    
    tdepth_k <- list()
    
    for (k in classes) {
      Yk <- reg_by_class[[as.character(k)]]
      if (is.null(Yk) || nrow(Yk) < 5) { tdepth_k[[as.character(k)]] <- NA_real_; next }
      
      # compute each curve's depth w.r.t. its own class sample
      dk <- rep(NA_real_, nrow(Yk))
      for (ii in seq_len(nrow(Yk))) {
        dk[ii] <- fm_depth_curve(Yk[ii, ], Yk, min_m_at_t = min_m_at_t)
      }
      
      tdepth_k[[as.character(k)]] <- get_typical(dk, typical_stat)
    }
    
    typical_depth[[region]] <- tdepth_k
    setTxtProgressBar(pb_typ, ri)
  }
  close(pb_typ)
  
  ## ----------------------------------------------------------
  ## Build normalized DD-features for TRAIN patients
  ## X_train: n_train x 3 matrix (normalized summed depths over regions)
  ## ----------------------------------------------------------
  train_patients <- sort(unique(train_df$RID))
  X_train <- matrix(NA_real_, nrow=length(train_patients), ncol=length(classes))
  colnames(X_train) <- c("d_CN","d_MCI","d_AD")
  y_train <- integer(length(train_patients))
  
  pb_tr <- txtProgressBar(min=0, max=length(train_patients), style=3)
  for (i in seq_along(train_patients)) {
    rid <- train_patients[i]
    subp <- train_df[train_df$RID == rid, ]
    subp <- subp[order(subp$month), ]
    y_train[i] <- subp$DIAGNOSIS[1]
    subp <- subp[1:min(nrow(subp), max_visits), ]
    
    S <- setNames(rep(0, length(classes)), as.character(classes))
    used <- FALSE
    
    for (region in region_names) {
      reg_models <- curves_train[[region]]
      reg_typ    <- typical_depth[[region]]
      if (is.null(reg_models) || is.null(reg_typ)) next
      
      yvec <- subp[[region]]
      if (all(is.na(yvec))) next
      idx_ok <- which(!is.na(yvec))
      if (length(idx_ok) < 2) next
      
      t <- subp$month[idx_ok]
      y <- yvec[idx_ok]
      
      if (use_pct_change) {
        if (!is.finite(y[1]) || y[1] == 0) next
        y <- 100 * (y - y[1]) / y[1]
      }
      
      y_grid <- interp_to_grid(t, y, grid)
      
      for (k in classes) {
        Y_train_k <- reg_models[[as.character(k)]]
        typ_k     <- reg_typ[[as.character(k)]]
        
        if (is.null(Y_train_k) || !is.finite(typ_k) || typ_k <= 0) next
        
        d_fk <- fm_depth_curve(y_grid, Y_train_k, min_m_at_t = min_m_at_t)
        if (is.finite(d_fk)) {
          # normalize depth by typical in-class self depth
          S[as.character(k)] <- S[as.character(k)] + d_fk / (typ_k + eps_typical)
          used <- TRUE
        }
      }
    }
    
    if (used) {
      X_train[i, ] <- c(S["1"], S["2"], S["4"])
    }
    setTxtProgressBar(pb_tr, i)
  }
  close(pb_tr)
  
  keep <- complete.cases(X_train) & is.finite(y_train)
  X_train <- X_train[keep, , drop=FALSE]
  y_train <- y_train[keep]
  
  ## ----------------------------------------------------------
  ## Predict TEST patients with normalized DD-features + kNN
  ## ----------------------------------------------------------
  test_patients <- sort(unique(test_df$RID))
  pb_te <- txtProgressBar(min=0, max=length(test_patients), style=3)
  
  fold_out <- data.frame(RID=integer(0), true=integer(0), pred=integer(0))
  majority_train <- as.numeric(names(sort(table(train_df$DIAGNOSIS), decreasing=TRUE))[1])
  
  for (i in seq_along(test_patients)) {
    rid <- test_patients[i]
    subp <- test_df[test_df$RID == rid, ]
    subp <- subp[order(subp$month), ]
    true <- subp$DIAGNOSIS[1]
    subp <- subp[1:min(nrow(subp), max_visits), ]
    
    S <- setNames(rep(0, length(classes)), as.character(classes))
    used <- FALSE
    
    for (region in region_names) {
      reg_models <- curves_train[[region]]
      reg_typ    <- typical_depth[[region]]
      if (is.null(reg_models) || is.null(reg_typ)) next
      
      yvec <- subp[[region]]
      if (all(is.na(yvec))) next
      idx_ok <- which(!is.na(yvec))
      if (length(idx_ok) < 2) next
      
      t <- subp$month[idx_ok]
      y <- yvec[idx_ok]
      
      if (use_pct_change) {
        if (!is.finite(y[1]) || y[1] == 0) next
        y <- 100 * (y - y[1]) / y[1]
      }
      
      y_grid <- interp_to_grid(t, y, grid)
      
      for (k in classes) {
        Y_train_k <- reg_models[[as.character(k)]]
        typ_k     <- reg_typ[[as.character(k)]]
        
        if (is.null(Y_train_k) || !is.finite(typ_k) || typ_k <= 0) next
        
        d_fk <- fm_depth_curve(y_grid, Y_train_k, min_m_at_t = min_m_at_t)
        if (is.finite(d_fk)) {
          S[as.character(k)] <- S[as.character(k)] + d_fk / (typ_k + eps_typical)
          used <- TRUE
        }
      }
    }
    
    if (!used) {
      pred <- majority_train
    } else {
      x_test <- matrix(c(S["1"], S["2"], S["4"]), nrow=1)
      colnames(x_test) <- colnames(X_train)
      pred <- knn(train = X_train, test = x_test, cl = y_train, k = k_nn)
      pred <- as.numeric(as.character(pred))
    }
    
    fold_out <- rbind(fold_out, data.frame(RID=rid, true=true, pred=pred))
    setTxtProgressBar(pb_te, i)
  }
  close(pb_te)
  
  all_preds[[fold]] <- fold_out
  
  cat("\nFold accuracy:", round(mean(fold_out$true == fold_out$pred), 4), "\n")
  cat("Confusion (rows=true, cols=pred):\n")
  print(table(True = label_class(fold_out$true), Pred = label_class(fold_out$pred)))
}

## ============================================================
## Overall performance
## ============================================================
preds <- do.call(rbind, all_preds)

cat("\n=============================\n")
cat("Overall CV accuracy:", round(mean(preds$true == preds$pred), 4), "\n")
cat("=============================\n")
cat("Overall confusion (rows=true, cols=pred):\n")
print(table(True = label_class(preds$true), Pred = label_class(preds$pred)))






library(dplyr)
library(class)

set.seed(123)

## -------------------------
## USER SETTINGS
## -------------------------
K_folds        <- 5
max_visits     <- 4
classes        <- c(1, 2, 4)

region_start   <- 17
region_cols    <- region_start:ncol(data)
# region_cols <- region_cols[1:60]   # optional speed-up for iteration

grid_step      <- 6
Tmax_use       <- 48       # cap horizon (important for stability)
min_m_at_t     <- 10
use_pct_change <- TRUE

k_nn           <- 30       # k in kNN
topN_regions   <- 40       # <<< try 20, 30, 40, 60, 80

## -------------------------
## BASIC CLEANUP
## -------------------------
data$DIAGNOSIS <- as.numeric(as.character(data$DIAGNOSIS))
data$month     <- as.numeric(as.character(data$month))

keep_cols <- c("RID","month","DIAGNOSIS", colnames(data)[region_cols])
df0 <- data[, keep_cols]
df0 <- df0[!is.na(df0$RID) & !is.na(df0$month) & !is.na(df0$DIAGNOSIS), ]
df0 <- df0[df0$DIAGNOSIS %in% classes, ]
colnames(df0) <- make.unique(colnames(df0))
region_names <- colnames(df0)[4:ncol(df0)]

label_class <- function(x) factor(x, levels=c(1,2,4), labels=c("CN","MCI","AD"))

## -------------------------
## Patient folds
## -------------------------
all_rids <- sort(unique(df0$RID))
fold_id  <- sample(rep(1:K_folds, length.out = length(all_rids)))
rid2fold <- data.frame(RID = all_rids, fold = fold_id)

## -------------------------
## Helper functions
## -------------------------
interp_to_grid <- function(t, y, grid) {
  if (length(t) < 2) return(rep(NA_real_, length(grid)))
  approx(x=t, y=y, xout=grid, method="linear", rule=1, ties="ordered")$y
}

fm_point_depth <- function(y0, y_train_vec) {
  y_train_vec <- y_train_vec[is.finite(y_train_vec)]
  m <- length(y_train_vec)
  if (!is.finite(y0) || m < 2) return(NA_real_)
  F <- mean(y_train_vec <= y0)
  1 - 2 * abs(F - 0.5)
}

fm_depth_curve <- function(y0_grid, Y_train_mat, min_m_at_t = 10) {
  n_grid <- length(y0_grid)
  d_t <- rep(NA_real_, n_grid)
  for (g in seq_len(n_grid)) {
    train_col <- Y_train_mat[, g]
    train_col <- train_col[is.finite(train_col)]
    if (length(train_col) < min_m_at_t) next
    d_t[g] <- fm_point_depth(y0_grid[g], train_col)
  }
  if (all(is.na(d_t))) return(NA_real_)
  mean(d_t, na.rm=TRUE)
}

## region separation score from class centroids in depth-space
sep_score <- function(mu1, mu2, mu4) {
  # average pairwise Euclidean distance between centroids
  d12 <- sqrt(sum((mu1 - mu2)^2))
  d14 <- sqrt(sum((mu1 - mu4)^2))
  d24 <- sqrt(sum((mu2 - mu4)^2))
  (d12 + d14 + d24) / 3
}

## ============================================================
## CROSS-VALIDATION
## ============================================================
all_preds <- list()

for (fold in 1:K_folds) {
  cat("\n=============================\n")
  cat("Fold", fold, "of", K_folds, "\n")
  cat("=============================\n")
  
  test_rids  <- rid2fold$RID[rid2fold$fold == fold]
  train_rids <- rid2fold$RID[rid2fold$fold != fold]
  
  train_df <- df0[df0$RID %in% train_rids, ]
  test_df  <- df0[df0$RID %in% test_rids, ]
  
  grid <- seq(0, Tmax_use, by = grid_step)
  
  ## ----------------------------------------------------------
  ## Build training curve matrices per region/class (<=4 visits)
  ## ----------------------------------------------------------
  curves_train <- setNames(vector("list", length(region_names)), region_names)
  
  pb_reg <- txtProgressBar(min=0, max=length(region_names), style=3)
  for (ri in seq_along(region_names)) {
    region <- region_names[ri]
    
    tmp <- train_df[, c("RID","month","DIAGNOSIS", region)]
    colnames(tmp)[4] <- "y"
    tmp <- tmp[!is.na(tmp$y), ]
    
    tmp <- tmp %>%
      arrange(RID, month) %>%
      group_by(RID) %>%
      slice_head(n = max_visits) %>%
      ungroup()
    
    reg_by_class <- list()
    for (k in classes) {
      subk <- tmp[tmp$DIAGNOSIS == k, ]
      rids_k <- sort(unique(subk$RID))
      if (length(rids_k) < 5) { reg_by_class[[as.character(k)]] <- NULL; next }
      
      Yk <- matrix(NA_real_, nrow=length(rids_k), ncol=length(grid))
      for (ii in seq_along(rids_k)) {
        rid <- rids_k[ii]
        si <- subk[subk$RID == rid, ]
        si <- si[order(si$month), ]
        t <- si$month; y <- si$y
        
        if (use_pct_change) {
          if (!is.finite(y[1]) || y[1] == 0) next
          y <- 100 * (y - y[1]) / y[1]
        }
        Yk[ii, ] <- interp_to_grid(t, y, grid)
      }
      reg_by_class[[as.character(k)]] <- Yk
    }
    
    curves_train[[region]] <- reg_by_class
    setTxtProgressBar(pb_reg, ri)
  }
  close(pb_reg)
  
  ## ----------------------------------------------------------
  ## Compute region weights from TRAIN data only
  ## Weight = separation between class centroids in depth-space
  ## ----------------------------------------------------------
  region_weight <- setNames(rep(NA_real_, length(region_names)), region_names)
  
  pb_w <- txtProgressBar(min=0, max=length(region_names), style=3)
  for (ri in seq_along(region_names)) {
    region <- region_names[ri]
    reg <- curves_train[[region]]
    if (is.null(reg)) { setTxtProgressBar(pb_w, ri); next }
    
    Y1 <- reg[["1"]]; Y2 <- reg[["2"]]; Y4 <- reg[["4"]]
    if (is.null(Y1) || is.null(Y2) || is.null(Y4)) { setTxtProgressBar(pb_w, ri); next }
    
    # compute mean depth vector per class for THIS region
    # For each curve, compute depths to each class, then average within true class
    
    mean_depth_vec <- function(Y_true, reg) {
      # returns (mean depth to CN, MCI, AD) for curves in Y_true
      out <- matrix(NA_real_, nrow=nrow(Y_true), ncol=3)
      for (ii in seq_len(nrow(Y_true))) {
        f <- Y_true[ii, ]
        out[ii, 1] <- fm_depth_curve(f, reg[["1"]], min_m_at_t)
        out[ii, 2] <- fm_depth_curve(f, reg[["2"]], min_m_at_t)
        out[ii, 3] <- fm_depth_curve(f, reg[["4"]], min_m_at_t)
      }
      colMeans(out, na.rm=TRUE)
    }
    
    mu1 <- mean_depth_vec(Y1, reg)
    mu2 <- mean_depth_vec(Y2, reg)
    mu4 <- mean_depth_vec(Y4, reg)
    
    if (any(!is.finite(c(mu1,mu2,mu4)))) {
      region_weight[region] <- NA_real_
    } else {
      region_weight[region] <- sep_score(mu1, mu2, mu4)
    }
    
    setTxtProgressBar(pb_w, ri)
  }
  close(pb_w)
  
  # keep topN regions by weight
  region_weight <- region_weight[is.finite(region_weight)]
  region_weight <- sort(region_weight, decreasing = TRUE)
  
  top_regions <- names(region_weight)[1:min(topN_regions, length(region_weight))]
  w_top <- region_weight[top_regions]
  
  cat("\nUsing top", length(top_regions), "regions by separation.\n")
  
  ## ----------------------------------------------------------
  ## Build DD-features (weighted sum over regions) for TRAIN
  ## ----------------------------------------------------------
  train_patients <- sort(unique(train_df$RID))
  X_train <- matrix(NA_real_, nrow=length(train_patients), ncol=3)
  colnames(X_train) <- c("d_CN","d_MCI","d_AD")
  y_train <- integer(length(train_patients))
  
  pb_tr <- txtProgressBar(min=0, max=length(train_patients), style=3)
  for (i in seq_along(train_patients)) {
    rid <- train_patients[i]
    subp <- train_df[train_df$RID == rid, ]
    subp <- subp[order(subp$month), ]
    y_train[i] <- subp$DIAGNOSIS[1]
    subp <- subp[1:min(nrow(subp), max_visits), ]
    
    S <- c(CN=0, MCI=0, AD=0)
    used <- FALSE
    
    for (region in top_regions) {
      reg <- curves_train[[region]]
      if (is.null(reg)) next
      
      yvec <- subp[[region]]
      if (all(is.na(yvec))) next
      idx_ok <- which(!is.na(yvec))
      if (length(idx_ok) < 2) next
      
      t <- subp$month[idx_ok]
      y <- yvec[idx_ok]
      
      if (use_pct_change) {
        if (!is.finite(y[1]) || y[1] == 0) next
        y <- 100 * (y - y[1]) / y[1]
      }
      y_grid <- interp_to_grid(t, y, grid)
      
      d1 <- fm_depth_curve(y_grid, reg[["1"]], min_m_at_t)
      d2 <- fm_depth_curve(y_grid, reg[["2"]], min_m_at_t)
      d4 <- fm_depth_curve(y_grid, reg[["4"]], min_m_at_t)
      
      if (all(!is.finite(c(d1,d2,d4)))) next
      
      w <- w_top[region]
      if (!is.finite(w)) next
      
      if (is.finite(d1)) S["CN"]  <- S["CN"]  + w * d1
      if (is.finite(d2)) S["MCI"] <- S["MCI"] + w * d2
      if (is.finite(d4)) S["AD"]  <- S["AD"]  + w * d4
      used <- TRUE
    }
    
    if (used) X_train[i, ] <- c(S["CN"], S["MCI"], S["AD"])
    setTxtProgressBar(pb_tr, i)
  }
  close(pb_tr)
  
  keep <- complete.cases(X_train) & is.finite(y_train)
  X_train <- X_train[keep, , drop=FALSE]
  y_train <- y_train[keep]
  
  ## ----------------------------------------------------------
  ## Predict TEST
  ## ----------------------------------------------------------
  test_patients <- sort(unique(test_df$RID))
  pb_te <- txtProgressBar(min=0, max=length(test_patients), style=3)
  
  fold_out <- data.frame(RID=integer(0), true=integer(0), pred=integer(0))
  majority_train <- as.numeric(names(sort(table(train_df$DIAGNOSIS), decreasing=TRUE))[1])
  
  for (i in seq_along(test_patients)) {
    rid <- test_patients[i]
    subp <- test_df[test_df$RID == rid, ]
    subp <- subp[order(subp$month), ]
    true <- subp$DIAGNOSIS[1]
    subp <- subp[1:min(nrow(subp), max_visits), ]
    
    S <- c(CN=0, MCI=0, AD=0)
    used <- FALSE
    
    for (region in top_regions) {
      reg <- curves_train[[region]]
      if (is.null(reg)) next
      
      yvec <- subp[[region]]
      if (all(is.na(yvec))) next
      idx_ok <- which(!is.na(yvec))
      if (length(idx_ok) < 2) next
      
      t <- subp$month[idx_ok]
      y <- yvec[idx_ok]
      
      if (use_pct_change) {
        if (!is.finite(y[1]) || y[1] == 0) next
        y <- 100 * (y - y[1]) / y[1]
      }
      y_grid <- interp_to_grid(t, y, grid)
      
      d1 <- fm_depth_curve(y_grid, reg[["1"]], min_m_at_t)
      d2 <- fm_depth_curve(y_grid, reg[["2"]], min_m_at_t)
      d4 <- fm_depth_curve(y_grid, reg[["4"]], min_m_at_t)
      
      if (all(!is.finite(c(d1,d2,d4)))) next
      
      w <- w_top[region]
      if (!is.finite(w)) next
      
      if (is.finite(d1)) S["CN"]  <- S["CN"]  + w * d1
      if (is.finite(d2)) S["MCI"] <- S["MCI"] + w * d2
      if (is.finite(d4)) S["AD"]  <- S["AD"]  + w * d4
      used <- TRUE
    }
    
    if (!used) {
      pred <- majority_train
    } else {
      x_test <- matrix(c(S["CN"], S["MCI"], S["AD"]), nrow=1)
      colnames(x_test) <- colnames(X_train)
      pred <- knn(train=X_train, test=x_test, cl=y_train, k=k_nn)
      pred <- as.numeric(as.character(pred))
    }
    
    fold_out <- rbind(fold_out, data.frame(RID=rid, true=true, pred=pred))
    setTxtProgressBar(pb_te, i)
  }
  close(pb_te)
  
  all_preds[[fold]] <- fold_out
  
  cat("\nFold accuracy:", round(mean(fold_out$true == fold_out$pred), 4), "\n")
  cat("Confusion (rows=true, cols=pred):\n")
  print(table(True = label_class(fold_out$true), Pred = label_class(fold_out$pred)))
}

preds <- do.call(rbind, all_preds)
cat("\n=============================\n")
cat("Overall CV accuracy:", round(mean(preds$true == preds$pred), 4), "\n")
cat("=============================\n")
cat("Overall confusion (rows=true, cols=pred):\n")
print(table(True = label_class(preds$true), Pred = label_class(preds$pred)))
