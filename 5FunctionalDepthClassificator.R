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
n_perm <- 1000
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


## Max-depth classifier ---------
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
#"Il modello mostra che la profondità funzionale è un potente fattore protettivo (OR < 0.001), indicando che la deviazione dalla centralità volumetrica è strettamente legata alla patologia. Tra i fattori di rischio, l'età incrementa le probabilità di AD dell'8.7% annuo, mentre il genotipo 4/4 rappresenta il rischio genetico maggiore, con un'odds 36 volte superiore rispetto al gruppo di riferimento 3/3."

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


# L'AUC di 0.83 indica un'ottima capacità discriminante. Significa che, prendendo a caso un paziente AD e un soggetto sano, il modello assegnerà una probabilità di malattia più alta al paziente AD nell'83% dei casi. È un valore molto solido per dati clinici complessi.


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


# Survival analysis -----------
library(dplyr)
library(survival)

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
