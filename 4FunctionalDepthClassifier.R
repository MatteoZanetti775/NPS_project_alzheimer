# Computation of Functional Depths -----------
library(dplyr)
library(ggplot2)
library(roahd)
library(tidyr)

datatot <- read.csv("data.csv")
data <- datatot[datatot$nav != 1, ]
data <- data %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) == 1) %>%
  ungroup() %>%
  mutate(DIAGNOSIS = case_when(
    DIAGNOSIS == 1 ~ "CN",
    DIAGNOSIS == 2 ~ "MCI",
    DIAGNOSIS == 4 ~ "AD",
    TRUE ~ as.character(DIAGNOSIS)
  ))

volume_names <- colnames(data)[17:133]
data[, volume_names] <- scale(data[, volume_names])
patients <- unique(data$RID)
grid <- seq(min(data$age, na.rm = TRUE), max(data$age, na.rm = TRUE), length.out = 50)

functional_list <- lapply(volume_names, function(vol) {
  t(sapply(patients, function(p) {
    p_data <- data[data$RID == p, ]
    p_data <- aggregate(p_data[[vol]] ~ p_data$age, FUN = mean)
    approx(p_data[[1]], p_data[[2]], xout = grid, rule = 2)$y
  }))
})

mf_data <- mfData(grid, functional_list)
depth_values <- multiMBD(mf_data)
depth_matrix <- data.frame(RID = patients, functional_depth = depth_values)
dim(depth_matrix)
head(depth_matrix)

depth_merged <- depth_matrix %>%
  left_join(data %>% distinct(RID, DIAGNOSIS), by = "RID")

dim(depth_merged)
head(depth_merged)

## Box-Plot --------------
ggplot(depth_merged, aes(x = DIAGNOSIS, y = functional_depth, fill = DIAGNOSIS)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Depth distribution (Stable patients)", y = "Functional Depth (MBD)")

## Permutation test -------
# CN vs AD
subset_CN_AD <- depth_merged %>% filter(DIAGNOSIS %in% c("CN", "AD"))
obs_diff <- mean(subset_CN_AD$functional_depth[subset_CN_AD$DIAGNOSIS == "CN"]) - 
  mean(subset_CN_AD$functional_depth[subset_CN_AD$DIAGNOSIS == "AD"])

perm_diffs <- replicate(n_perm, {
  shuffled_labels <- sample(subset_CN_AD$DIAGNOSIS)
  mean(subset_CN_AD$functional_depth[shuffled_labels == "CN"]) - 
    mean(subset_CN_AD$functional_depth[shuffled_labels == "AD"])
})
p_value <- mean(abs(perm_diffs) >= abs(obs_diff))
cat("Permutation test p-value (CN vs AD stable):", p_value, "\n")

# CN vs MCI
subset_CN_MCI <- depth_merged %>% filter(DIAGNOSIS %in% c("CN", "MCI"))
obs_diff <- mean(subset_CN_MCI$functional_depth[subset_CN_MCI$DIAGNOSIS == "CN"]) - 
  mean(subset_CN_MCI$functional_depth[subset_CN_MCI$DIAGNOSIS == "MCI"])

perm_diffs <- replicate(n_perm, {
  shuffled_labels <- sample(subset_CN_MCI$DIAGNOSIS)
  mean(subset_CN_MCI$functional_depth[shuffled_labels == "CN"]) - 
    mean(subset_CN_MCI$functional_depth[shuffled_labels == "MCI"])
})
p_value <- mean(abs(perm_diffs) >= abs(obs_diff))
cat("Permutation test p-value (CN vs MCI stable):", p_value, "\n")

# MCI vs AD
subset_MCI_AD <- depth_merged %>% filter(DIAGNOSIS %in% c("MCI", "AD"))
obs_diff <- mean(subset_MCI_AD$functional_depth[subset_MCI_AD$DIAGNOSIS == "MCI"]) - 
  mean(subset_MCI_AD$functional_depth[subset_MCI_AD$DIAGNOSIS == "AD"])

perm_diffs <- replicate(n_perm, {
  shuffled_labels <- sample(subset_MCI_AD$DIAGNOSIS)
  mean(subset_MCI_AD$functional_depth[shuffled_labels == "MCI"]) - 
    mean(subset_MCI_AD$functional_depth[shuffled_labels == "AD"])
})
p_value <- mean(abs(perm_diffs) >= abs(obs_diff))
cat("Permutation test p-value (MCI vs AD stable):", p_value, "\n")


#write.csv(depth_merged, "depth_stable_patients.csv", row.names = FALSE)

depth_merged = read.csv("depth_stable_patients.csv")

# Classification -----------

## Max-depth classsifier --------
library(dplyr)
library(roahd)
library(tidyr)

datatot <- read.csv("data.csv")
volume_names <- colnames(datatot)[17:133]

data_stable <- datatot %>%
  filter(nav != 1) %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) == 1) %>%
  ungroup() %>%
  mutate(DIAGNOSIS = case_when(DIAGNOSIS == 1 ~ "CN", DIAGNOSIS == 2 ~ "MCI", DIAGNOSIS == 4 ~ "AD", TRUE ~ as.character(DIAGNOSIS)))

data_conv <- datatot %>%
  filter(nav != 1) %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) != 1) %>%
  ungroup() %>%
  mutate(DIAGNOSIS_FINAL = case_when(last(DIAGNOSIS) == 1 ~ "CN", last(DIAGNOSIS) == 2 ~ "MCI", last(DIAGNOSIS) == 4 ~ "AD", TRUE ~ as.character(last(DIAGNOSIS))))

means_stable <- colMeans(data_stable[, volume_names], na.rm = TRUE)
sds_stable <- apply(data_stable[, volume_names], 2, sd, na.rm = TRUE)

data_stable[, volume_names] <- scale(data_stable[, volume_names], center = means_stable, scale = sds_stable)
data_conv[, volume_names] <- scale(data_conv[, volume_names], center = means_stable, scale = sds_stable)

grid <- seq(min(datatot$age, na.rm = TRUE), max(datatot$age, na.rm = TRUE), length.out = 50)

create_mf_data <- function(df, p_list, v_names, g_points) {
  lapply(v_names, function(vol) {
    t(sapply(p_list, function(p) {
      p_data <- df[df$RID == p, ]
      p_data <- aggregate(p_data[[vol]] ~ p_data$age, FUN = mean)
      approx(p_data[[1]], p_data[[2]], xout = g_points, rule = 2)$y
    }))
  })
}

patients_stable <- unique(data_stable$RID)
patients_conv <- unique(data_conv$RID)

mf_stable <- mfData(grid, create_mf_data(data_stable, patients_stable, volume_names, grid))
mf_conv <- mfData(grid, create_mf_data(data_conv, patients_conv, volume_names, grid))

mbd_custom <- function(X, P) {
  if(is.vector(X)) X <- matrix(X, nrow = 1)
  n <- nrow(P); L <- ncol(P); depths <- numeric(nrow(X))
  for (i in 1:nrow(X)) {
    count <- 0
    for (t in 1:L) {
      n_a <- sum(P[, t] <= X[i, t]); n_b <- sum(P[, t] >= X[i, t])
      count <- count + (n_a * n_b)
    }
    depths[i] <- count / (L * n * (n + 1) / 2)
  }
  return(depths)
}

calc_multi_depth_final <- function(mfd_to_test, mfd_reference) {
  depth_matrix <- matrix(0, nrow = mfd_to_test$N, ncol = mfd_to_test$L)
  for(i in 1:mfd_to_test$L) {
    depth_matrix[, i] <- mbd_custom(mfd_to_test$fDList[[i]]$values, mfd_reference$fDList[[i]]$values)
  }
  rowMeans(depth_matrix)
}

mf_cn <- mf_stable[which(patients_stable %in% (data_stable %>% filter(DIAGNOSIS == "CN") %>% pull(RID)))]
mf_mci <- mf_stable[which(patients_stable %in% (data_stable %>% filter(DIAGNOSIS == "MCI") %>% pull(RID)))]
mf_ad <- mf_stable[which(patients_stable %in% (data_stable %>% filter(DIAGNOSIS == "AD") %>% pull(RID)))]

# STABLE ACCURACY
class_stable <- data.frame(
  Actual = (data_stable %>% distinct(RID, DIAGNOSIS) %>% arrange(match(RID, patients_stable)) %>% pull(DIAGNOSIS)),
  d_CN = calc_multi_depth_final(mf_stable, mf_cn),
  d_MCI = calc_multi_depth_final(mf_stable, mf_mci),
  d_AD = calc_multi_depth_final(mf_stable, mf_ad)
) %>% mutate(Pred = case_when(d_CN >= d_MCI & d_CN >= d_AD ~ "CN", d_MCI >= d_CN & d_MCI >= d_AD ~ "MCI", TRUE ~ "AD"))

# CONVERTERS ACCURACY
class_conv <- data.frame(
  Actual = (data_conv %>% group_by(RID) %>% summarise(DIAG_FINAL = first(DIAGNOSIS_FINAL)) %>% arrange(match(RID, patients_conv)) %>% pull(DIAG_FINAL)),
  d_CN = calc_multi_depth_final(mf_conv, mf_cn),
  d_MCI = calc_multi_depth_final(mf_conv, mf_mci),
  d_AD = calc_multi_depth_final(mf_conv, mf_ad)
) %>% mutate(Pred = case_when(d_CN >= d_MCI & d_CN >= d_AD ~ "CN", d_MCI >= d_CN & d_MCI >= d_AD ~ "MCI", TRUE ~ "AD"))

m_stable <- table(class_stable$Actual, class_stable$Pred)
m_conv <- table(class_conv$Actual, class_conv$Pred)

cat("ACCURACY STABLE PATIENTS (Training):\n")
print(m_stable)
cat("Accuracy:", sum(diag(m_stable))/sum(m_stable), "\n\n")

cat("ACCURACY CONVERTERS (Test - Final Diagnosis Prediction):\n")
print(m_conv)
cat("Accuracy:", sum(diag(m_conv))/sum(m_conv), "\n")


## Max-depth classifier: relevant volumes only --------
library(dplyr)
library(roahd)
library(tidyr)

datatot <- read.csv("data.csv")

volume_names <- c(
  "Left.Hippocampus.Volume..ST29SV.", "Right.Hippocampus.Volume..ST88SV.",
  "Left.Entorhinal.Volume..ST24CV.", "Right.Entorhinal.Volume..ST83CV.",
  "Left.Amygdala.Volume..ST12SV.", "Right.Amygdala.Volume..ST76SV.",
  "Left.Parahippocampal.Volume..ST44CV.", "Right.Parahippocampal.Volume..ST103CV.",
  "Left.Fusiform.Volume..ST26CV.", "Right.Fusiform.Volume..ST85CV.",
  "Left.Lateral.Ventricle.Volume..ST37SV.", "Right.Lateral.Ventricle.Volume..ST96SV.",
  "Left.Posterior.Cingulate.Volume..ST50CV.", "Right.Posterior.Cingulate.Volume..ST109CV.",
  "Left.Precuneus.Volume..ST54CV.", "Right.Precuneus.Volume..ST111CV.",
  "Left.Middle.Temporal.Volume..ST40CV.", "Right.Middle.Temporal.Volume..ST99CV."
)

data_stable <- datatot %>%
  filter(nav != 1) %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) == 1) %>%
  ungroup() %>%
  mutate(DIAGNOSIS = case_when(
    DIAGNOSIS == 1 ~ "CN", 
    DIAGNOSIS == 2 ~ "MCI", 
    DIAGNOSIS == 4 ~ "AD", 
    TRUE ~ as.character(DIAGNOSIS)
  ))

data_conv <- datatot %>%
  filter(nav != 1) %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) != 1) %>%
  ungroup() %>%
  mutate(DIAGNOSIS_FINAL = case_when(
    last(DIAGNOSIS) == 1 ~ "CN", 
    last(DIAGNOSIS) == 2 ~ "MCI", 
    last(DIAGNOSIS) == 4 ~ "AD", 
    TRUE ~ as.character(last(DIAGNOSIS))
  ))

means_stable <- colMeans(data_stable[, volume_names], na.rm = TRUE)
sds_stable <- apply(data_stable[, volume_names], 2, sd, na.rm = TRUE)

data_stable[, volume_names] <- scale(data_stable[, volume_names], center = means_stable, scale = sds_stable)
data_conv[, volume_names] <- scale(data_conv[, volume_names], center = means_stable, scale = sds_stable)

grid <- seq(min(datatot$age, na.rm = TRUE), max(datatot$age, na.rm = TRUE), length.out = 50)

create_mf_data <- function(df, p_list, v_names, g_points) {
  lapply(v_names, function(vol) {
    t(sapply(p_list, function(p) {
      p_data <- df[df$RID == p, ]
      p_data <- aggregate(p_data[[vol]] ~ p_data$age, FUN = mean)
      approx(p_data[[1]], p_data[[2]], xout = g_points, rule = 2)$y
    }))
  })
}

patients_stable <- unique(data_stable$RID)
patients_conv <- unique(data_conv$RID)

mf_stable <- mfData(grid, create_mf_data(data_stable, patients_stable, volume_names, grid))
mf_conv <- mfData(grid, create_mf_data(data_conv, patients_conv, volume_names, grid))

mbd_custom <- function(X, P) {
  if(is.vector(X)) X <- matrix(X, nrow = 1)
  n <- nrow(P)
  L <- ncol(P)
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
  return(depths)
}

calc_multi_depth_final <- function(mfd_to_test, mfd_reference) {
  depth_matrix <- matrix(0, nrow = mfd_to_test$N, ncol = mfd_to_test$L)
  for(i in 1:mfd_to_test$L) {
    depth_matrix[, i] <- mbd_custom(mfd_to_test$fDList[[i]]$values, mfd_reference$fDList[[i]]$values)
  }
  rowMeans(depth_matrix)
}

mf_cn <- mf_stable[which(patients_stable %in% (data_stable %>% filter(DIAGNOSIS == "CN") %>% pull(RID)))]
mf_mci <- mf_stable[which(patients_stable %in% (data_stable %>% filter(DIAGNOSIS == "MCI") %>% pull(RID)))]
mf_ad <- mf_stable[which(patients_stable %in% (data_stable %>% filter(DIAGNOSIS == "AD") %>% pull(RID)))]

class_stable <- data.frame(
  Actual = (data_stable %>% distinct(RID, DIAGNOSIS) %>% arrange(match(RID, patients_stable)) %>% pull(DIAGNOSIS)),
  d_CN = calc_multi_depth_final(mf_stable, mf_cn),
  d_MCI = calc_multi_depth_final(mf_stable, mf_mci),
  d_AD = calc_multi_depth_final(mf_stable, mf_ad)
) %>% mutate(Pred = case_when(d_CN >= d_MCI & d_CN >= d_AD ~ "CN", d_MCI >= d_CN & d_MCI >= d_AD ~ "MCI", TRUE ~ "AD"))

class_conv <- data.frame(
  Actual = (data_conv %>% group_by(RID) %>% summarise(DIAG_FINAL = first(DIAGNOSIS_FINAL)) %>% arrange(match(RID, patients_conv)) %>% pull(DIAG_FINAL)),
  d_CN = calc_multi_depth_final(mf_conv, mf_cn),
  d_MCI = calc_multi_depth_final(mf_conv, mf_mci),
  d_AD = calc_multi_depth_final(mf_conv, mf_ad)
) %>% mutate(Pred = case_when(d_CN >= d_MCI & d_CN >= d_AD ~ "CN", d_MCI >= d_CN & d_MCI >= d_AD ~ "MCI", TRUE ~ "AD"))

m_stable <- table(class_stable$Actual, class_stable$Pred)
m_conv <- table(class_conv$Actual, class_conv$Pred)

cat("ACCURACY STABLE PATIENTS (Training):\n")
print(m_stable)
cat("Accuracy:", sum(diag(m_stable))/sum(m_stable), "\n\n")

cat("ACCURACY CONVERTERS (Test - Final Diagnosis Prediction):\n")
print(m_conv)
cat("Accuracy:", sum(diag(m_conv))/sum(m_conv), "\n")

## Max-depth classifier: 2 groups ---------
library(dplyr)
library(roahd)
library(tidyr)

datatot <- read.csv("data.csv")
volume_names <- colnames(datatot)[17:133]

data_stable <- datatot %>%
  filter(nav != 1) %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) == 1) %>%
  ungroup() %>%
  mutate(DIAGNOSIS = case_when(DIAGNOSIS == 1 ~ "CN", DIAGNOSIS == 4 ~ "AD", TRUE ~ "OTHER")) %>%
  filter(DIAGNOSIS %in% c("CN", "AD"))

data_conv <- datatot %>%
  filter(nav != 1) %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) != 1) %>%
  ungroup() %>%
  mutate(DIAGNOSIS_FINAL = case_when(last(DIAGNOSIS) == 1 ~ "CN", last(DIAGNOSIS) == 4 ~ "AD", TRUE ~ "OTHER")) %>%
  filter(DIAGNOSIS_FINAL %in% c("CN", "AD"))

means_stable <- colMeans(data_stable[, volume_names], na.rm = TRUE)
sds_stable <- apply(data_stable[, volume_names], 2, sd, na.rm = TRUE)
data_stable[, volume_names] <- scale(data_stable[, volume_names], center = means_stable, scale = sds_stable)
data_conv[, volume_names] <- scale(data_conv[, volume_names], center = means_stable, scale = sds_stable)

grid <- seq(min(datatot$age, na.rm = TRUE), max(datatot$age, na.rm = TRUE), length.out = 50)

create_mf_data <- function(df, p_list, v_names, g_points) {
  lapply(v_names, function(vol) {
    t(sapply(p_list, function(p) {
      p_data <- df[df$RID == p, ]
      p_data <- aggregate(p_data[[vol]] ~ p_data$age, FUN = mean)
      approx(p_data[[1]], p_data[[2]], xout = g_points, rule = 2)$y
    }))
  })
}

patients_stable <- unique(data_stable$RID)
patients_conv <- unique(data_conv$RID)

mf_stable <- mfData(grid, create_mf_data(data_stable, patients_stable, volume_names, grid))
mf_conv <- mfData(grid, create_mf_data(data_conv, patients_conv, volume_names, grid))

mbd_custom <- function(X, P) {
  if(is.vector(X)) X <- matrix(X, nrow = 1)
  n <- nrow(P); L <- ncol(P); depths <- numeric(nrow(X))
  for (i in 1:nrow(X)) {
    count <- 0
    for (t in 1:L) {
      n_a <- sum(P[, t] <= X[i, t]); n_b <- sum(P[, t] >= X[i, t])
      count <- count + (n_a * n_b)
    }
    depths[i] <- count / (L * n * (n + 1) / 2)
  }
  return(depths)
}

calc_multi_depth_final <- function(mfd_to_test, mfd_reference) {
  depth_matrix <- matrix(0, nrow = mfd_to_test$N, ncol = mfd_to_test$L)
  for(i in 1:mfd_to_test$L) {
    depth_matrix[, i] <- mbd_custom(mfd_to_test$fDList[[i]]$values, mfd_reference$fDList[[i]]$values)
  }
  rowMeans(depth_matrix)
}

mf_cn <- mf_stable[which(patients_stable %in% (data_stable %>% filter(DIAGNOSIS == "CN") %>% pull(RID)))]
mf_ad <- mf_stable[which(patients_stable %in% (data_stable %>% filter(DIAGNOSIS == "AD") %>% pull(RID)))]

class_stable <- data.frame(
  Actual = (data_stable %>% distinct(RID, DIAGNOSIS) %>% arrange(match(RID, patients_stable)) %>% pull(DIAGNOSIS)),
  d_CN = calc_multi_depth_final(mf_stable, mf_cn),
  d_AD = calc_multi_depth_final(mf_stable, mf_ad)
) %>% mutate(Pred = if_else(d_CN >= d_AD, "CN", "AD"))

class_conv <- data.frame(
  Actual = (data_conv %>% group_by(RID) %>% summarise(F = first(DIAGNOSIS_FINAL)) %>% arrange(match(RID, patients_conv)) %>% pull(F)),
  d_CN = calc_multi_depth_final(mf_conv, mf_cn),
  d_AD = calc_multi_depth_final(mf_conv, mf_ad)
) %>% mutate(Pred = if_else(d_CN >= d_AD, "CN", "AD"))

m_stable <- table(class_stable$Actual, class_stable$Pred)
m_conv <- table(class_conv$Actual, class_conv$Pred)

cat("BINARY ACCURACY STABLE:\n")
print(m_stable)
cat("Accuracy:", sum(diag(m_stable))/sum(m_stable), "\n\n")

cat("BINARY ACCURACY CONVERTERS:\n")
print(m_conv)
cat("Accuracy:", sum(diag(m_conv))/sum(m_conv), "\n")



## Logistic regression: depth + covariates ------------
data_reg <- depth_merged %>%
  left_join(datatot %>%
              group_by(RID) %>%
              summarise(age = min(age, na.rm = TRUE),
                        GENOTYPE = first(GENOTYPE),
                        .groups = "drop"),
            by = "RID") %>%
  filter(DIAGNOSIS %in% c("CN", "AD")) %>%
  mutate(Outcome = if_else(DIAGNOSIS == "AD", 1, 0))

data_reg$GENOTYPE <- factor(data_reg$GENOTYPE)
data_reg$GENOTYPE <- relevel(data_reg$GENOTYPE, ref = "3/3")

model <- glm(Outcome ~ functional_depth + age + factor(GENOTYPE),
             data = data_reg,
             family = binomial)

summary(model)
table(data_reg$GENOTYPE, data_reg$Outcome)
obs_coeff <- summary(model)$coefficients["functional_depth", "Estimate"]
n_perm <- 1000
set.seed(123)

perm_coeffs <- replicate(n_perm, {
  data_reg$shuffled_outcome <- sample(data_reg$Outcome)
  perm_mod <- glm(shuffled_outcome ~ functional_depth + age + factor(GENOTYPE),
                  data = data_reg, family = binomial)
  if("functional_depth" %in% rownames(summary(perm_mod)$coefficients)) {
    summary(perm_mod)$coefficients["functional_depth", "Estimate"]
  } else {
    NA
  }
})

perm_coeffs <- perm_coeffs[!is.na(perm_coeffs)]
p_val_coeff <- mean(abs(perm_coeffs) >= abs(obs_coeff))

cat("P-value permutazionale per l'effetto della Profondità:", p_val_coeff, "\n")

### odds ratio --------
exp(coef(model))

### Accuracy --------
data_reg$prob_pred <- predict(model, type = "response")
data_reg$class_pred <- if_else(data_reg$prob_pred > 0.5, 1, 0)

conf_matrix <- table(Actual = data_reg$Outcome, Predicted = data_reg$class_pred)
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)

cat("MATRICE DI CONFUSIONE (Logistica):\n")
print(conf_matrix)
cat("\nAccuratezza complessiva:", round(accuracy, 3), "\n")

library(pROC)

roc_obj <- roc(data_reg$Outcome, data_reg$prob_pred)
cat("AUC del modello:", round(auc(roc_obj), 3), "\n")

plot(roc_obj, 
     main = paste("Curva ROC (CN vs AD) - AUC:", round(auc(roc_obj), 3)), 
     col = "#2c7bb6", 
     lwd = 3)
abline(a = 0, b = 1, lty = 2, col = "grey")

#### Accuracy on data_conv -------
data_test_meta <- data_conv %>%
  group_by(RID) %>%
  summarise(
    age = min(age, na.rm = TRUE),
    GENOTYPE = first(GENOTYPE),
    Outcome = if_else(last(DIAGNOSIS) == 4, 1, 0),
    .groups = "drop"
  )

mf_test <- mfData(grid, create_mf_data(data_conv, data_test_meta$RID, volume_names, grid))

data_test_meta$functional_depth <- calc_multi_depth_final(mf_test, mf_cn_ref)

data_test_meta <- data_test_meta %>% filter(GENOTYPE %in% unique(data_reg_all$GENOTYPE))
data_test_meta$GENOTYPE <- factor(data_test_meta$GENOTYPE, levels = levels(data_reg_all$GENOTYPE))

data_test_meta$prob_pred <- predict(model, newdata = data_test_meta, type = "response")
data_test_meta$class_pred <- if_else(data_test_meta$prob_pred > 0.5, 1, 0)

conf_matrix_test <- table(Actual = data_test_meta$Outcome, Predicted = data_test_meta$class_pred)
accuracy_test <- sum(diag(conf_matrix_test)) / sum(conf_matrix_test)

print(conf_matrix_test)
cat("Accuracy on Test Set (Converters):", round(accuracy_test, 3), "\n")
library(pROC)

roc_conv <- roc(data_test_meta$Outcome, data_test_meta$prob_pred)
auc_conv <- auc(roc_conv)

cat("AUC sui Converters (modello allenato su tutti gli stabili):", round(auc_conv, 3), "\n")

plot(roc_conv, 
     main = paste("ROC Converters - AUC:", round(auc_conv, 3)), 
     col = "#d95f02", 
     lwd = 3)
abline(a = 0, b = 1, lty = 2, col = "grey")

### Computing Indice di Youden ------------
coords_opt <- coords(roc_obj, "best", ret = "threshold", best.method = "youden")
best_threshold <- coords_opt$threshold

cat("The optimal threshold is:", round(best_threshold, 3), "\n")

data_reg$class_pred_youden <- if_else(data_reg$prob_pred > best_threshold, 1, 0)

conf_matrix_youden <- table(Actual = data_reg$Outcome, Predicted = data_reg$class_pred_youden)
accuracy_youden <- sum(diag(conf_matrix_youden)) / sum(conf_matrix_youden)
data_test_meta$class_pred_youden <- if_else(data_test_meta$prob_pred > best_threshold, 1, 0)
acc_conv_youden <- mean(data_test_meta$class_pred_youden == data_test_meta$Outcome)

print(table(Actual = data_test_meta$Outcome, Predicted = data_test_meta$class_pred_youden))
cat("Accuracy on Converters with Youden threshold:", round(acc_conv_youden, 3), "\n")

cat("MATRICE DI CONFUSIONE (Youden Threshold):\n")
print(conf_matrix_youden)
cat("\nAccuratezza ottimizzata:", round(accuracy_youden, 3), "\n")

# Logistic regression: depth + covariates, test set (stable patients) -------
library(dplyr)
library(roahd)
library(tidyr)
library(pROC)

datatot <- read.csv("data.csv")
volume_names <- colnames(datatot)[17:133]

data_stable_all <- datatot %>%
  filter(nav != 1) %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) == 1) %>%
  ungroup() %>%
  mutate(DIAGNOSIS = case_when(
    DIAGNOSIS == 1 ~ "CN",
    DIAGNOSIS == 2 ~ "MCI",
    DIAGNOSIS == 4 ~ "AD",
    TRUE ~ as.character(DIAGNOSIS)
  ))

data_conv <- datatot %>%
  filter(nav != 1) %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) != 1) %>%
  ungroup() %>%
  mutate(DIAGNOSIS_FINAL = if_else(last(DIAGNOSIS) == 4, "AD", "CN"))

means_stable <- colMeans(data_stable_all[, volume_names], na.rm = TRUE)
sds_stable <- apply(data_stable_all[, volume_names], 2, sd, na.rm = TRUE)

data_stable_all[, volume_names] <- scale(data_stable_all[, volume_names], center = means_stable, scale = sds_stable)
data_conv[, volume_names] <- scale(data_conv[, volume_names], center = means_stable, scale = sds_stable)

grid <- seq(min(datatot$age, na.rm = TRUE), max(datatot$age, na.rm = TRUE), length.out = 50)

create_mf_data <- function(df, p_list, v_names, g_points) {
  lapply(v_names, function(vol) {
    t(sapply(p_list, function(p) {
      p_data <- df[df$RID == p, ]
      p_data <- aggregate(p_data[[vol]] ~ p_data$age, FUN = mean)
      approx(p_data[[1]], p_data[[2]], xout = g_points, rule = 2)$y
    }))
  })
}

patients_stable_all <- unique(data_stable_all$RID)
mf_stable_all <- mfData(grid, create_mf_data(data_stable_all, patients_stable_all, volume_names, grid))
mf_cn_ref <- mf_stable_all[which(patients_stable_all %in% (data_stable_all %>% filter(DIAGNOSIS == "CN") %>% pull(RID)))]

mbd_custom <- function(X, P) {
  if(is.vector(X)) X <- matrix(X, nrow = 1)
  n <- nrow(P); L <- ncol(P); depths <- numeric(nrow(X))
  for (i in 1:nrow(X)) {
    count <- 0
    for (t in 1:L) {
      n_a <- sum(P[, t] <= X[i, t]); n_b <- sum(P[, t] >= X[i, t])
      count <- count + (n_a * n_b)
    }
    depths[i] <- count / (L * n * (n + 1) / 2)
  }
  return(depths)
}

calc_multi_depth_final <- function(mfd_to_test, mfd_reference) {
  depth_matrix <- matrix(0, nrow = mfd_to_test$N, ncol = mfd_to_test$L)
  for(i in 1:mfd_to_test$L) {
    depth_matrix[, i] <- mbd_custom(mfd_to_test$fDList[[i]]$values, mfd_reference$fDList[[i]]$values)
  }
  rowMeans(depth_matrix)
}

depth_all_stable <- data.frame(
  RID = patients_stable_all,
  functional_depth = calc_multi_depth_final(mf_stable_all, mf_cn_ref)
)

data_reg_full <- depth_all_stable %>%
  left_join(data_stable_all %>% distinct(RID, DIAGNOSIS), by = "RID") %>%
  left_join(datatot %>% group_by(RID) %>% summarise(age = min(age), GENOTYPE = first(GENOTYPE), .groups = "drop"), by = "RID") %>%
  filter(DIAGNOSIS %in% c("CN", "AD")) %>%
  mutate(Outcome = if_else(DIAGNOSIS == "AD", 1, 0))

data_reg_full$GENOTYPE <- factor(data_reg_full$GENOTYPE)
data_reg_full$GENOTYPE <- relevel(data_reg_full$GENOTYPE, ref = "3/3")

set.seed(123)
train_idx <- sample(1:nrow(data_reg_full), size = 0.8 * nrow(data_reg_full))
data_train <- data_reg_full[train_idx, ]
data_stable_test <- data_reg_full[-train_idx, ]

#Remove levels not present in training
levels_in_train <- unique(data_train$GENOTYPE)
data_stable_test <- data_stable_test %>% filter(GENOTYPE %in% levels_in_train)

model <- glm(Outcome ~ functional_depth + age + GENOTYPE, data = data_train, family = binomial)

data_stable_test$prob_pred <- predict(model, newdata = data_stable_test, type = "response")
accuracy_stable_test <- mean((data_stable_test$prob_pred > 0.5) == data_stable_test$Outcome)

data_test_meta <- data_conv %>%
  group_by(RID) %>%
  summarise(age = min(age, na.rm = TRUE), GENOTYPE = first(GENOTYPE), Outcome = if_else(last(DIAGNOSIS) == 4, 1, 0), .groups = "drop")

mf_conv_test <- mfData(grid, create_mf_data(data_conv, data_test_meta$RID, volume_names, grid))
data_test_meta$functional_depth <- calc_multi_depth_final(mf_conv_test, mf_cn_ref)

data_test_meta <- data_test_meta %>% filter(GENOTYPE %in% levels_in_train)
data_test_meta$GENOTYPE <- factor(data_test_meta$GENOTYPE, levels = levels(data_train$GENOTYPE))

data_test_meta$prob_pred <- predict(model, newdata = data_test_meta, type = "response")
accuracy_conv_test <- mean((data_test_meta$prob_pred > 0.5) == data_test_meta$Outcome)

cat("ACCURACY STABLE TEST (Stable 20%):", round(accuracy_stable_test, 3), "\n")
cat("ACCURACY CONVERTERS TEST:", round(accuracy_conv_test, 3), "\n")

conf_stable <- table(Actual = data_stable_test$Outcome, Predicted = as.numeric(data_stable_test$prob_pred > 0.5))
conf_conv <- table(Actual = data_test_meta$Outcome, Predicted = as.numeric(data_test_meta$prob_pred > 0.5))

print("Confusion Matrix - Stable Test:")
print(conf_stable)
print("Confusion Matrix - Converters Test:")
print(conf_conv)


roc_stable <- roc(data_stable_test$Outcome, data_stable_test$prob_pred)
auc_stable <- auc(roc_stable)

roc_conv <- roc(data_test_meta$Outcome, data_test_meta$prob_pred)
auc_conv <- auc(roc_conv)

cat("AUC Stable Test Set:", round(auc_stable, 3), "\n")
cat("AUC Converters Test Set:", round(auc_conv, 3), "\n")

par(mfrow=c(1,2))

plot(roc_stable, 
     main = paste("ROC Stable Test\nAUC:", round(auc_stable, 3)), 
     col = "#2c7bb6", lwd = 3)
abline(a = 0, b = 1, lty = 2, col = "grey")

plot(roc_conv, 
     main = paste("ROC Converters Test\nAUC:", round(auc_conv, 3)), 
     col = "#d95f02", lwd = 3)
abline(a = 0, b = 1, lty = 2, col = "grey")

par(mfrow=c(1,1))


### Computing Indice di Youden ------------
coords_opt_conv <- coords(roc_conv, "best", ret = "threshold", best.method = "youden")
best_threshold_conv <- coords_opt_conv$threshold

cat("The optimal threshold for Converters is:", round(best_threshold_conv, 3), "\n")

data_stable_test$class_pred_youden <- if_else(data_stable_test$prob_pred > best_threshold_conv, 1, 0)
data_test_meta$class_pred_youden <- if_else(data_test_meta$prob_pred > best_threshold_conv, 1, 0)

conf_stable_youden <- table(Actual = data_stable_test$Outcome, Predicted = data_stable_test$class_pred_youden)
conf_conv_youden <- table(Actual = data_test_meta$Outcome, Predicted = data_test_meta$class_pred_youden)

acc_stable_youden <- sum(diag(conf_stable_youden)) / sum(conf_stable_youden)
acc_conv_youden <- sum(diag(conf_conv_youden)) / sum(conf_conv_youden)

cat("\n--- FINAL RESULTS WITH YOUDEN THRESHOLD ---\n")
cat("Accuracy Stable Test (Youden):", round(acc_stable_youden, 3), "\n")
cat("Accuracy Converters Test (Youden):", round(acc_conv_youden, 3), "\n")

print("Confusion Matrix - Stable Test (Youden):")
print(conf_stable_youden)
print("Confusion Matrix - Converters Test (Youden):")
print(conf_conv_youden)


# Logistic regression IP(CN, MCI, AD) ~ functional_depth + covariates --------
data_reg_all <- depth_merged %>%
  left_join(datatot %>%
              group_by(RID) %>%
              summarise(age = min(age, na.rm = TRUE),
                        GENOTYPE = first(GENOTYPE),
                        .groups = "drop"),
            by = "RID") %>%
  filter(DIAGNOSIS %in% c("CN", "MCI", "AD")) %>%
  mutate(Outcome = if_else(DIAGNOSIS == "AD", 1, 0))

data_reg_all$GENOTYPE <- factor(data_reg_all$GENOTYPE)
data_reg_all$GENOTYPE <- relevel(data_reg_all$GENOTYPE, ref = "3/3")

model_all <- glm(Outcome ~ functional_depth + age + GENOTYPE,
                 data = data_reg_all,
                 family = binomial)

summary(model_all)

data_reg_all$prob_pred <- predict(model_all, type = "response")
data_reg_all$class_pred <- if_else(data_reg_all$prob_pred > 0.5, 1, 0)

conf_matrix_all <- table(Actual = data_reg_all$Outcome, Predicted = data_reg_all$class_pred)
accuracy_all <- sum(diag(conf_matrix_all)) / sum(conf_matrix_all)

print(conf_matrix_all)
cat("Accuracy with MCI included:", round(accuracy_all, 3), "\n")

roc_all <- roc(data_reg_all$Outcome, data_reg_all$prob_pred)
cat("AUC with MCI included:", round(auc(roc_all), 3), "\n")

plot(roc_all, main = "ROC Curve: (CN+MCI) vs AD", col = "#d95f02", lwd = 3) #buona capacità discriminante



### Accuracy on data_conv --------
data_test_meta <- data_conv %>%
  group_by(RID) %>%
  summarise(
    age = min(age, na.rm = TRUE),
    GENOTYPE = first(GENOTYPE),
    Outcome = if_else(last(DIAGNOSIS) == 4, 1, 0),
    .groups = "drop"
  )

mf_test <- mfData(grid, create_mf_data(data_conv, data_test_meta$RID, volume_names, grid))
data_test_meta$functional_depth <- calc_multi_depth_final(mf_test, mf_cn_ref)

data_test_meta <- data_test_meta %>% filter(GENOTYPE %in% unique(data_reg_all$GENOTYPE))
data_test_meta$GENOTYPE <- factor(data_test_meta$GENOTYPE, levels = levels(data_reg_all$GENOTYPE))

data_test_meta$prob_pred <- predict(model_all, newdata = data_test_meta, type = "response")
data_test_meta$class_pred <- if_else(data_test_meta$prob_pred > 0.5, 1, 0)
data_test_meta$class_pred <- if_else(data_test_meta$prob_pred > 0.35, 1, 0)

conf_matrix_test <- table(Actual = data_test_meta$Outcome, Predicted = data_test_meta$class_pred)
accuracy_test <- sum(diag(conf_matrix_test)) / sum(conf_matrix_test)

print(conf_matrix_test)
cat("Accuracy on Test Set (Converters):", round(accuracy_test, 3), "\n")
library(pROC)
data_test_meta$prob_pred_all <- predict(model_all, newdata = data_test_meta, type = "response")
roc_conv_all <- roc(data_test_meta$Outcome, data_test_meta$prob_pred_all)
auc_conv_all <- auc(roc_conv_all)
cat("AUC sui Converters (Modello AD vs non-AD):", round(auc_conv_all, 3), "\n")
plot(roc_conv_all, 
     main = paste("ROC Converters (AD vs non-AD) - AUC:", round(auc_conv_all, 3)), 
     col = "#d95f02", 
     lwd = 3)
abline(a = 0, b = 1, lty = 2, col = "grey")

### Computing Indice di Youden ------------
library(pROC)

roc_conv <- roc(data_test_meta$Outcome, data_test_meta$prob_pred, quiet = TRUE)
auc_conv <- auc(roc_conv)

coords_opt <- coords(roc_conv, "best", ret = "threshold", best.method = "youden")
best_threshold <- coords_opt$threshold

cat("The optimal Youden threshold is:", round(best_threshold, 3), "\n")

data_test_meta$class_pred_youden <- if_else(data_test_meta$prob_pred > best_threshold, 1, 0)

conf_matrix_youden <- table(Actual = data_test_meta$Outcome, Predicted = data_test_meta$class_pred_youden)
accuracy_youden <- sum(diag(conf_matrix_youden)) / sum(conf_matrix_youden)

print(conf_matrix_youden)
cat("Accuracy on Test Set (Youden):", round(accuracy_youden, 3), "\n")
cat("AUC on Test Set:", round(auc_conv, 3), "\n")

plot(roc_conv, 
     main = paste("ROC Test Set (Converters) - AUC:", round(auc_conv, 3)), 
     col = "#d95f02", 
     lwd = 3)
abline(a = 0, b = 1, lty = 2, col = "grey")


# Logistic regression: IP(CN, MCI, AD) ~ functional_depth + covariates, stable patients (test) --------
library(dplyr)
library(roahd)
library(tidyr)
library(pROC)

datatot <- read.csv("data.csv")
volume_names <- colnames(datatot)[17:133]

data_stable_all <- datatot %>%
  filter(nav != 1) %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) == 1) %>%
  ungroup() %>%
  mutate(DIAGNOSIS = case_when(
    DIAGNOSIS == 1 ~ "CN",
    DIAGNOSIS == 2 ~ "MCI",
    DIAGNOSIS == 4 ~ "AD",
    TRUE ~ as.character(DIAGNOSIS)
  ))

data_conv <- datatot %>%
  filter(nav != 1) %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) != 1) %>%
  ungroup()

means_stable <- colMeans(data_stable_all[, volume_names], na.rm = TRUE)
sds_stable <- apply(data_stable_all[, volume_names], 2, sd, na.rm = TRUE)

data_stable_all[, volume_names] <- scale(data_stable_all[, volume_names], center = means_stable, scale = sds_stable)
data_conv[, volume_names] <- scale(data_conv[, volume_names], center = means_stable, scale = sds_stable)

grid <- seq(min(datatot$age, na.rm = TRUE), max(datatot$age, na.rm = TRUE), length.out = 50)

create_mf_data <- function(df, p_list, v_names, g_points) {
  lapply(v_names, function(vol) {
    t(sapply(p_list, function(p) {
      p_data <- df[df$RID == p, ]
      p_data <- aggregate(p_data[[vol]] ~ p_data$age, FUN = mean)
      approx(p_data[[1]], p_data[[2]], xout = g_points, rule = 2)$y
    }))
  })
}

patients_stable_all <- unique(data_stable_all$RID)
mf_stable_all <- mfData(grid, create_mf_data(data_stable_all, patients_stable_all, volume_names, grid))
mf_cn_ref <- mf_stable_all[which(patients_stable_all %in% (data_stable_all %>% filter(DIAGNOSIS == "CN") %>% pull(RID)))]

mbd_custom <- function(X, P) {
  if(is.vector(X)) X <- matrix(X, nrow = 1)
  n <- nrow(P); L <- ncol(P); depths <- numeric(nrow(X))
  for (i in 1:nrow(X)) {
    count <- 0
    for (t in 1:L) {
      n_a <- sum(P[, t] <= X[i, t]); n_b <- sum(P[, t] >= X[i, t])
      count <- count + (n_a * n_b)
    }
    depths[i] <- count / (L * n * (n + 1) / 2)
  }
  return(depths)
}

calc_multi_depth_final <- function(mfd_to_test, mfd_reference) {
  depth_matrix <- matrix(0, nrow = mfd_to_test$N, ncol = mfd_to_test$L)
  for(i in 1:mfd_to_test$L) {
    depth_matrix[, i] <- mbd_custom(mfd_to_test$fDList[[i]]$values, mfd_reference$fDList[[i]]$values)
  }
  rowMeans(depth_matrix)
}

depth_all_stable <- data.frame(
  RID = patients_stable_all,
  functional_depth = calc_multi_depth_final(mf_stable_all, mf_cn_ref)
)

data_reg_all <- depth_all_stable %>%
  left_join(data_stable_all %>% distinct(RID, DIAGNOSIS), by = "RID") %>%
  left_join(datatot %>% group_by(RID) %>% summarise(age = min(age), GENOTYPE = first(GENOTYPE), .groups = "drop"), by = "RID") %>%
  filter(DIAGNOSIS %in% c("CN", "MCI", "AD")) %>%
  mutate(Outcome = if_else(DIAGNOSIS == "AD", 1, 0))

data_reg_all$GENOTYPE <- factor(data_reg_all$GENOTYPE)
data_reg_all$GENOTYPE <- relevel(data_reg_all$GENOTYPE, ref = "3/3")

set.seed(123)
train_idx <- sample(1:nrow(data_reg_all), size = 0.8 * nrow(data_reg_all))
data_train <- data_reg_all[train_idx, ]
data_stable_test <- data_reg_all[-train_idx, ]

levels_in_train <- unique(data_train$GENOTYPE)
data_stable_test <- data_stable_test %>% filter(GENOTYPE %in% levels_in_train)

model_all <- glm(Outcome ~ functional_depth + age + GENOTYPE, data = data_train, family = binomial)

data_stable_test$prob_pred <- predict(model_all, newdata = data_stable_test, type = "response")
acc_stable <- mean((data_stable_test$prob_pred > 0.5) == data_stable_test$Outcome)
roc_stable <- roc(data_stable_test$Outcome, data_stable_test$prob_pred, quiet = TRUE)

data_test_meta <- data_conv %>%
  group_by(RID) %>%
  summarise(age = min(age, na.rm = TRUE), GENOTYPE = first(GENOTYPE), Outcome = if_else(last(DIAGNOSIS) == 4, 1, 0), .groups = "drop")

mf_conv_test <- mfData(grid, create_mf_data(data_conv, data_test_meta$RID, volume_names, grid))
data_test_meta$functional_depth <- calc_multi_depth_final(mf_conv_test, mf_cn_ref)

data_test_meta <- data_test_meta %>% filter(GENOTYPE %in% levels_in_train)
data_test_meta$GENOTYPE <- factor(data_test_meta$GENOTYPE, levels = levels(data_train$GENOTYPE))

data_test_meta$prob_pred <- predict(model_all, newdata = data_test_meta, type = "response")
acc_conv <- mean((data_test_meta$prob_pred > 0.5) == data_test_meta$Outcome)

cat("--- RESULTS ---\n")
cat("Accuracy Stable Test:", round(acc_stable, 3), "\n")
cat("AUC Stable Test:", round(auc(roc_stable), 3), "\n")
cat("Accuracy Converters Test:", round(acc_conv, 3), "\n")

if(length(unique(data_test_meta$Outcome)) > 1) {
  roc_conv <- roc(data_test_meta$Outcome, data_test_meta$prob_pred, quiet = TRUE)
  cat("AUC Converters Test:", round(auc(roc_conv), 3), "\n")
  
  par(mfrow=c(1,2))
  plot(roc_stable, main = paste("ROC Stable\nAUC:", round(auc(roc_stable), 3)), col = "#2c7bb6", lwd = 3)
  plot(roc_conv, main = paste("ROC Converters\nAUC:", round(auc(roc_conv), 3)), col = "#d95f02", lwd = 3)
  par(mfrow=c(1,1))
} else {
  cat("AUC Converters: N/A (Only 1 class present)\n")
  plot(roc_stable, main = paste("ROC Stable\nAUC:", round(auc(roc_stable), 3)), col = "#2c7bb6", lwd = 3)
}


### Computing Indice di Youden ------------
library(pROC)

roc_conv <- roc(data_test_meta$Outcome, data_test_meta$prob_pred, quiet = TRUE)
coords_opt <- coords(roc_conv, "best", ret = "threshold", best.method = "youden")
best_threshold <- coords_opt$threshold

cat("The optimal Youden threshold identified is:", round(best_threshold, 3), "\n\n")

data_stable_test$class_pred_youden <- if_else(data_stable_test$prob_pred > best_threshold, 1, 0)
conf_stable_youden <- table(Actual = data_stable_test$Outcome, Predicted = data_stable_test$class_pred_youden)
acc_stable_youden <- sum(diag(conf_stable_youden)) / sum(conf_stable_youden)

data_test_meta$class_pred_youden <- if_else(data_test_meta$prob_pred > best_threshold, 1, 0)
conf_conv_youden <- table(Actual = data_test_meta$Outcome, Predicted = data_test_meta$class_pred_youden)
acc_conv_youden <- sum(diag(conf_conv_youden)) / sum(conf_conv_youden)

cat("--- FINAL RESULTS (Threshold:", round(best_threshold, 3), ") ---\n")

cat("\n[STABLE TEST SET]\n")
print(conf_stable_youden)
cat("Accuracy Stable (Youden):", round(acc_stable_youden, 3), "\n")
cat("AUC Stable:", round(auc(roc_stable), 3), "\n")

cat("\n[CONVERTERS TEST SET]\n")
print(conf_conv_youden)
cat("Accuracy Converters (Youden):", round(acc_conv_youden, 3), "\n")
cat("AUC Converters:", round(auc(roc_conv), 3), "\n")

par(mfrow=c(1,2))
plot(roc_stable, main = paste("Stable Test - AUC:", round(auc(roc_stable), 3)), col = "#2c7bb6", lwd = 3)
abline(h = coords(roc_stable, x = best_threshold, input = "threshold", ret = "sensitivity"), col = "red", lty = 3)

plot(roc_conv, main = paste("Conv Test - AUC:", round(auc(roc_conv), 3)), col = "#d95f02", lwd = 3)
abline(h = coords(roc_conv, x = best_threshold, input = "threshold", ret = "sensitivity"), col = "red", lty = 3)
par(mfrow=c(1,1))



# Multinomial regression  ----------
library(dplyr)
library(nnet)
library(tidyr)
library(pROC)

data_multinom <- depth_merged %>%
  left_join(datatot %>%
              group_by(RID) %>%
              summarise(age = min(age, na.rm = TRUE),
                        GENOTYPE = first(GENOTYPE),
                        .groups = "drop"),
            by = "RID") %>%
  filter(DIAGNOSIS %in% c("CN", "MCI", "AD"))

data_multinom$DIAGNOSIS <- factor(data_multinom$DIAGNOSIS)
data_multinom$DIAGNOSIS <- relevel(data_multinom$DIAGNOSIS, ref = "CN")
data_multinom$GENOTYPE <- factor(data_multinom$GENOTYPE)
data_multinom$GENOTYPE <- relevel(data_multinom$GENOTYPE, ref = "3/3")

model_multi <- multinom(DIAGNOSIS ~ functional_depth + age + GENOTYPE, 
                        data = data_multinom)

summary(model_multi)

probs_train <- predict(model_multi, type = "probs")
pred_train <- predict(model_multi, type = "class")

conf_matrix_multi <- table(Actual = data_multinom$DIAGNOSIS, Predicted = pred_train)
accuracy_multi <- sum(diag(conf_matrix_multi)) / sum(conf_matrix_multi)

cat("MULTINOMIAL CONFUSION MATRIX:\n")
print(conf_matrix_multi)
cat("\nAccuracy:", round(accuracy_multi, 3), "\n")

exp(coef(model_multi))

data_test_conv <- data_conv %>%
  group_by(RID) %>%
  summarise(
    age = min(age, na.rm = TRUE),
    GENOTYPE = first(GENOTYPE),
    Outcome_Final = last(DIAGNOSIS_FINAL),
    .groups = "drop"
  )


# correction
patients_stable <- unique(data_stable$RID)

mf_stable <- mfData(
  grid,
  create_mf_data(data_stable, patients_stable, volume_names, grid)
)

mf_cn_ref <- mf_stable[
  which(patients_stable %in%
          (data_stable %>%
             filter(DIAGNOSIS == "CN") %>%
             pull(RID)))
]


mf_test <- mfData(grid, create_mf_data(data_conv, data_test_conv$RID, volume_names, grid))
data_test_conv$functional_depth <- calc_multi_depth_final(mf_test, mf_cn_ref)

data_test_conv$GENOTYPE <- factor(data_test_conv$GENOTYPE, levels = levels(data_multinom$GENOTYPE))
data_test_conv$Outcome_Final <- factor(data_test_conv$Outcome_Final, levels = levels(data_multinom$DIAGNOSIS))

probs_test <- predict(model_multi, newdata = data_test_conv, type = "probs")
pred_test <- predict(model_multi, newdata = data_test_conv, type = "class")

conf_matrix_test <- table(Actual = data_test_conv$Outcome_Final, Predicted = pred_test)
accuracy_test <- sum(diag(conf_matrix_test)) / sum(conf_matrix_test)

cat("\nACCURACY ON CONVERTERS (Multinomial):\n")
print(conf_matrix_test)
cat("Accuracy:", round(accuracy_test, 3), "\n")


# logistic regression: other covariates ----------
library(dplyr)
library(nnet)
library(pROC)
library(tidyr)
library(roahd)

datatot <- read.csv("data.csv")
volume_names <- colnames(datatot)[17:133]
grid <- seq(min(datatot$age, na.rm = TRUE), max(datatot$age, na.rm = TRUE), length.out = 50)

create_mf_data <- function(df, p_list, v_names, g_points) {
  lapply(v_names, function(vol) {
    t(sapply(p_list, function(p) {
      p_data <- df[df$RID == p, ]
      p_data <- aggregate(p_data[[vol]] ~ p_data$age, FUN = mean)
      approx(p_data[[1]], p_data[[2]], xout = g_points, rule = 2)$y
    }))
  })
}

mbd_custom <- function(X, P) {
  if(is.vector(X)) X <- matrix(X, nrow = 1)
  n <- nrow(P); L <- ncol(P); depths <- numeric(nrow(X))
  for (i in 1:nrow(X)) {
    count <- 0
    for (t in 1:L) {
      n_a <- sum(P[, t] <= X[i, t]); n_b <- sum(P[, t] >= X[i, t])
      count <- count + (n_a * n_b)
    }
    depths[i] <- count / (L * n * (n + 1) / 2)
  }
  return(depths)
}

calc_multi_depth_final <- function(mfd_to_test, mfd_reference) {
  depth_matrix <- matrix(0, nrow = mfd_to_test$N, ncol = mfd_to_test$L)
  for(i in 1:mfd_to_test$L) {
    depth_matrix[, i] <- mbd_custom(mfd_to_test$fDList[[i]]$values, mfd_reference$fDList[[i]]$values)
  }
  rowMeans(depth_matrix)
}

data_stable_all <- datatot %>%
  filter(nav != 1) %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) == 1) %>%
  ungroup() %>%
  mutate(DIAGNOSIS = case_when(
    DIAGNOSIS == 1 ~ "CN",
    DIAGNOSIS == 2 ~ "MCI",
    DIAGNOSIS == 4 ~ "AD",
    TRUE ~ as.character(DIAGNOSIS)
  ))

data_conv <- datatot %>%
  filter(nav != 1) %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) != 1) %>%
  ungroup()

means_stable <- colMeans(data_stable_all[, volume_names], na.rm = TRUE)
sds_stable <- apply(data_stable_all[, volume_names], 2, sd, na.rm = TRUE)
data_stable_all[, volume_names] <- scale(data_stable_all[, volume_names], center = means_stable, scale = sds_stable)
data_conv[, volume_names] <- scale(data_conv[, volume_names], center = means_stable, scale = sds_stable)

patients_stable_all <- unique(data_stable_all$RID)
mf_stable_all <- mfData(grid, create_mf_data(data_stable_all, patients_stable_all, volume_names, grid))
mf_cn_ref <- mf_stable_all[which(patients_stable_all %in% (data_stable_all %>% filter(DIAGNOSIS == "CN") %>% pull(RID)))]

depth_all_stable <- data.frame(
  RID = patients_stable_all,
  functional_depth = calc_multi_depth_final(mf_stable_all, mf_cn_ref)
)

metadata_covariates <- datatot %>%
  group_by(RID) %>%
  summarise(age_baseline = min(age, na.rm = TRUE),
            GENOTYPE = first(GENOTYPE),
            DXDEP = first(DXDEP),
            .groups = "drop") %>%
  mutate(GENOTYPE = factor(GENOTYPE),
         GENOTYPE = relevel(GENOTYPE, ref = "3/3"),
         DXDEP = factor(DXDEP),
         age_scaled = as.numeric(scale(age_baseline)))

data_train_stable <- depth_all_stable %>%
  inner_join(data_stable_all %>% distinct(RID, DIAGNOSIS), by = "RID") %>%
  inner_join(metadata_covariates, by = "RID")

rid_conv <- unique(data_conv$RID)
mf_conv_test <- mfData(grid, create_mf_data(data_conv, rid_conv, volume_names, grid))
depth_conv_val <- calc_multi_depth_final(mf_conv_test, mf_cn_ref)

data_test_converters <- data.frame(RID = rid_conv, functional_depth = depth_conv_val) %>%
  left_join(metadata_covariates, by = "RID") %>%
  left_join(data_conv %>%
              group_by(RID) %>%
              summarise(Outcome_Final = if_else(last(DIAGNOSIS) == 4, "AD", "CN"), .groups = "drop"),
            by = "RID")

cat("\n--- BINARY REGRESSION (Stable CN vs AD) ---\n")
data_bin_train <- data_train_stable %>% filter(DIAGNOSIS %in% c("CN", "AD")) %>%
  mutate(Outcome = if_else(DIAGNOSIS == "AD", 1, 0))

model_bin <- glm(Outcome ~ functional_depth + age_scaled + GENOTYPE + DXDEP,
                 data = data_bin_train, family = binomial)

data_test_converters$prob_bin <- predict(model_bin, newdata = data_test_converters, type = "response")
roc_bin <- roc(if_else(data_test_converters$Outcome_Final == "AD", 1, 0), data_test_converters$prob_bin, quiet = TRUE)
bin_acc <- mean((data_test_converters$prob_bin > 0.5) == (data_test_converters$Outcome_Final == "AD"))

print(summary(model_bin))
cat("AUC on Converters:", round(auc(roc_bin), 3), "\n")
cat("Accuracy on Converters (Binary):", round(bin_acc, 3), "\n")

cat("\n--- MULTINOMIAL REGRESSION (Stable CN vs MCI vs AD) ---\n")
data_multi_train <- data_train_stable %>% filter(DIAGNOSIS %in% c("CN", "MCI", "AD"))
data_multi_train$DIAGNOSIS <- factor(data_multi_train$DIAGNOSIS, levels = c("CN", "MCI", "AD"))

model_multi <- multinom(DIAGNOSIS ~ functional_depth + age_scaled + GENOTYPE + DXDEP,
                        data = data_multi_train, trace = FALSE)

data_test_converters$pred_multi <- predict(model_multi, newdata = data_test_converters)
conf_matrix_multi <- table(Actual = data_test_converters$Outcome_Final, Predicted = data_test_converters$pred_multi)
multi_acc <- sum(diag(conf_matrix_multi)) / sum(conf_matrix_multi)

print(conf_matrix_multi)
cat("Accuracy on Converters (Multinomial):", round(multi_acc, 3), "\n")
print(exp(coef(model_multi)))


# (Survival analysis) -----------
library(dplyr)
library(survival)
datatot = read.csv("data.csv")
data_surv_prep <- datatot %>%
  filter(!is.na(DIAGNOSIS)) %>%
  group_by(RID) %>%
  summarise(
    converter = if_else(any(DIAGNOSIS == 4), 1, 0),
    time = if_else(any(DIAGNOSIS == 4), 
                   min(month[DIAGNOSIS == 4], na.rm = TRUE), 
                   max(month, na.rm = TRUE)),
    age_baseline = min(age, na.rm = TRUE),
    GENOTYPE = first(GENOTYPE),
    .groups = "drop"
  )

data_survival <- data_surv_prep %>%
  left_join(depth_matrix, by = "RID") %>%
  filter(!is.na(functional_depth))

data_survival$GENOTYPE <- factor(data_survival$GENOTYPE)
data_survival$GENOTYPE <- relevel(data_survival$GENOTYPE, ref = "3/3")

cox_model <- coxph(Surv(time, converter) ~ functional_depth + age_baseline + GENOTYPE, 
                   data = data_survival)

summary(cox_model)

hazard_ratios <- exp(coef(cox_model))
print("Hazard Ratios:")
print(hazard_ratios)
