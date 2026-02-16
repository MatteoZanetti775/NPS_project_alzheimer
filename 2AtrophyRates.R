library(lme4)
library(quantreg)
library(lqmm)
library(dplyr)
library(ggplot2)
library(tidyr)
library(brms)
library(caret)
# 0) newdata0 Quantile regression ---------------
datatot = read.csv("data.csv")
idno = which(datatot$nav !=1)
data = datatot[idno,]
data = ungroup(filter(group_by(data, RID), n_distinct(DIAGNOSIS) == 1))
data$month = data$month/12
colnames(data) [5] = "years"
volumes = colnames(data)[17:133]
id1 = which(data$nav ==1)
u_rids = unique(data$RID)
n = length(u_rids)
newdata = data.frame(RID = u_rids)
newdata$DIAGNOSIS <- NA
newdata$GENOTYPE <- NA
newdata[volumes] <- NA
for (i in 1:n){
  curr_id = u_rids[i]
  patient_data = data[data$RID == curr_id, ]
  newdata[i, "DIAGNOSIS"] = patient_data$DIAGNOSIS[1]
  newdata[i, "GENOTYPE"] = patient_data$GENOTYPE[1]
  for (name in volumes){
    f = as.formula(paste(name, "~ years"))
    qlm = rq(f, data = patient_data, tau=0.5) 
    slope = coef(qlm)["years"]          # Tasso assoluto (mm3/anno)
    intercept = coef(qlm)["(Intercept)"] # Volume basale stimato (mm3)
    if(intercept != 0) {
      pct_change = (slope / intercept)*100
      newdata[i, name] = pct_change
    } else {
      newdata[i, name] = NA
    }
  }
  if(i %% 50 == 0) print(paste("Processati", i, "pazienti"))
}
newdata0 = newdata
View(newdata)
colnames(newdata)











# Quantile Regression changing --------------
datatot = read.csv("data.csv")
idno = which(datatot$nav !=1)
data = datatot[idno,]
data = ungroup(filter(group_by(data, RID), n_distinct(DIAGNOSIS) != 1))
data2 <- filter(
  group_by(data, RID),
  n_distinct(DIAGNOSIS) == 2 & first(DIAGNOSIS) != last(DIAGNOSIS)
)
data2 <- mutate(
  data2,
  DIAGNOSIS_suppl = paste0(first(DIAGNOSIS), "_to_", last(DIAGNOSIS))
)
data2 <- relocate(data2, DIAGNOSIS_suppl, .after = DIAGNOSIS)
# Seleziona solo i cambi desiderati: 1→2, 2→1, 2→4
data2 <- filter(
  group_by(data2, RID),
  DIAGNOSIS_suppl %in% c("1_to_2", "2_to_1", "2_to_4")
)
# Filtra pazienti con esattamente 3 diagnosi distinte
data3 <- filter(
  group_by(data, RID),
  n_distinct(DIAGNOSIS) == 3
)
data3 <- mutate(
  data3,
  DIAGNOSIS_suppl = paste(DIAGNOSIS, collapse = "_to_")
)
data3 <- filter(
  group_by(data3, RID),
  DIAGNOSIS_suppl == "1_to_2_to_4"
)
data3 <- relocate(data3, DIAGNOSIS_suppl, .after = DIAGNOSIS)
data <- bind_rows(data2, data3)

data$month = data$month/12
colnames(data) [5] = "years"
volumes = colnames(data)[18:134]
u_rids = unique(data$RID)
n = length(u_rids)
newdata = data.frame(RID = u_rids)
newdata$DIAGNOSIS <- NA
newdata$GENOTYPE <- NA
newdata[volumes] <- NA
for (i in 1:n){
  curr_id = u_rids[i]
  patient_data = data[data$RID == curr_id, ]
  newdata[i, "DIAGNOSIS"] = patient_data$DIAGNOSIS_suppl[1]
  newdata[i, "GENOTYPE"] = patient_data$GENOTYPE[1]
  for (name in volumes){
    f = as.formula(paste(name, "~ years"))
    qlm = rq(f, data = patient_data, tau=0.5) 
    slope = coef(qlm)["years"]          #Tasso assoluto (mm3/anno)
    intercept = coef(qlm)["(Intercept)"] #Volume basale stimato (mm3)
    if(intercept != 0) {
      pct_change = (slope / intercept)*100
      newdata[i, name] = pct_change
    } else {
      newdata[i, name] = NA
    }
  }
  if(i %% 50 == 0) print(paste("Processed", i, "patients"))
}
View(newdata)
colnames(newdata)
















# 1) newdata1 LMM ------------------------
datatot = read.csv("data.csv")
idno = which(datatot$nav != 1)
data = datatot[idno,]
data = ungroup(filter(group_by(data, RID), n_distinct(DIAGNOSIS) == 1))
data$month = data$month / 12
colnames(data)[5] = "years"

volumes = colnames(data)[17:133]
u_rids = unique(data$RID)

info = data %>%
  group_by(RID) %>%
  summarise(
    DIAGNOSIS = first(DIAGNOSIS),
    GENOTYPE = first(GENOTYPE)
  )

newdata = data.frame(RID = u_rids)
newdata = left_join(newdata, info, by = "RID")
newdata[volumes] <- NA

for (name in volumes) {
  
  f = as.formula(paste(name, "~ years + (years || RID)"))
  lmm = lmer(f, data = data, REML = TRUE)
  
  re = ranef(lmm)$RID
  
  slopes = fixef(lmm)["years"] + re[, "years"]
  intercepts = fixef(lmm)["(Intercept)"] + re[, "(Intercept)"]
  apc = (slopes / intercepts) * 100
  
  rid_names <- rownames(re)
  idx <- match(rid_names, newdata$RID)
  ok <- !is.na(idx)
  
  newdata[idx[ok], name] <- apc[ok]
}
newdata1 = newdata











# 2) newdata2 LMM: genotype -----------
datatot = read.csv("data.csv")
idno = which(datatot$nav != 1)
data = datatot[idno,]
data = data %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) == 1) %>%
  ungroup()

data$month = data$month / 12
colnames(data)[5] = "years"

data$DIAGNOSIS <- as.factor(data$DIAGNOSIS)
data$GENOTYPE <- as.factor(data$GENOTYPE)

volumes = colnames(data)[17:133]
u_rids = unique(data$RID)

info = data %>%
  group_by(RID) %>%
  summarise(
    DIAGNOSIS = first(DIAGNOSIS),
    GENOTYPE = first(GENOTYPE)
  )

newdata = data.frame(RID = u_rids)
newdata = left_join(newdata, info, by = "RID")
newdata[volumes] <- NA

for (name in volumes) {
  
  f = as.formula(paste(name, "~ years * GENOTYPE + (years | RID)"))
  
  tryCatch({
    model = lmer(f, data = data, REML = TRUE)
    
    pred_frame <- info
    
    pred_frame$years <- 0
    y0 <- predict(model, newdata = pred_frame)
    
    pred_frame$years <- 1
    y1 <- predict(model, newdata = pred_frame)
    
    slopes <- y1 - y0
    intercepts <- y0
    
    apc_values <- (slopes / intercepts) * 100
    
    idx <- match(pred_frame$RID, newdata$RID)
    newdata[idx, name] <- apc_values
    
  }, error = function(e) { 
    cat("Errore nel volume:", name, "\n") 
  })
}
newdata2 = newdata






# 3) newdata3 LMM: + genotype + depression + nav + volTot-------------
datatot = read.csv("data.csv")
idno = which(datatot$nav != 1)
data = datatot[idno,]
data = data %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) == 1) %>%
  ungroup()

data$month = data$month / 12
colnames(data)[5] = "years"

data$DIAGNOSIS <- as.factor(data$DIAGNOSIS)
data$GENOTYPE <- as.factor(data$GENOTYPE)
data$DXDEP <- as.factor(data$DXDEP)

volumes = colnames(data)[17:133]
u_rids = unique(data$RID)

col_supra <- "Supratentorial.Total.Volume..Aseg...ST128SV."
hs_vector <- numeric(length(u_rids))
for (i in 1:length(u_rids)) {
  sub_p <- data[data$RID == u_rids[i], ]
  sub_p <- sub_p[order(sub_p$years), ]
  valid_vals <- sub_p[[col_supra]][!is.na(sub_p[[col_supra]])]
  hs_vector[i] <- if(length(valid_vals) > 0) valid_vals[1] else NA
}
hs_df <- data.frame(RID = u_rids, HeadSize = hs_vector)
data <- merge(data, hs_df, by = "RID", all.x = TRUE)

info = data %>%
  group_by(RID) %>%
  summarise(
    DIAGNOSIS = first(DIAGNOSIS),
    GENOTYPE = first(GENOTYPE),
    HeadSize = first(HeadSize),
    nav = first(nav)
  )

newdata = data.frame(RID = u_rids)
newdata = left_join(newdata, info, by = "RID")
newdata[volumes] <- NA

for (name in volumes) {
  f <- as.formula(paste(name, "~ years * GENOTYPE + DXDEP + HeadSize + nav + (years | RID)"))
  
  tryCatch({
    model = lmer(f, data = data, REML = TRUE)
    
    pred_frame <- info
    pred_frame$DXDEP <- factor(0, levels = levels(data$DXDEP))
    
    pred_frame$years <- 0
    y0 <- predict(model, newdata = pred_frame)
    
    pred_frame$years <- 1
    y1 <- predict(model, newdata = pred_frame)
    
    slopes <- y1 - y0
    intercepts <- y0
    
    apc_values <- (slopes / intercepts) * 100
    
    idx <- match(pred_frame$RID, newdata$RID)
    newdata[idx, name] <- apc_values
    
  }, error = function(e) { 
    cat("Errore nel volume:", name, "\n") 
  })
}
newdata = newdata[, !(names(newdata) %in% c("nav", "HeadSize"))]
newdata3 = newdata











# 4) newdata4 Bayesian: genotype + depression + nav + VolTot <-- using "4Tests" ---------
datatot = read.csv("data.csv")
idno = which(datatot$nav != 1)
data = datatot[idno,]
data = data %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) == 1) %>%
  ungroup()

data$month = data$month / 12
colnames(data)[5] = "years"

data$DIAGNOSIS <- as.factor(data$DIAGNOSIS)
data$GENOTYPE <- as.factor(data$GENOTYPE)
data$DXDEP <- as.factor(data$DXDEP)

volumes = colnames(data)[17:133]
volumes = volumes[c(40, 47, 55, 45, 98)] #40 55 47 98 45 based on cleandata3
u_rids = unique(data$RID)

col_supra <- "Supratentorial.Total.Volume..Aseg...ST128SV."
hs_vector <- numeric(length(u_rids))
for (i in 1:length(u_rids)) {
  sub_p <- data[data$RID == u_rids[i], ]
  sub_p <- sub_p[order(sub_p$years), ]
  valid_vals <- sub_p[[col_supra]][!is.na(sub_p[[col_supra]])]
  hs_vector[i] <- if(length(valid_vals) > 0) valid_vals[1] else NA
}
hs_df <- data.frame(RID = u_rids, HeadSize = hs_vector)
data <- merge(data, hs_df, by = "RID", all.x = TRUE)

info = data %>%
  group_by(RID) %>%
  summarise(
    DIAGNOSIS = first(DIAGNOSIS),
    GENOTYPE = first(GENOTYPE),
    HeadSize = first(HeadSize),
    nav = first(nav)
  )

newdata = data.frame(RID = u_rids)
newdata = left_join(newdata, info, by = "RID")
newdata[volumes] <- NA

options(mc.cores = 2) 

cat("--- Avvio compilazione modello base ---\n")
first_vol <- volumes[1]
formula_base <- as.formula(paste(first_vol, "~ years*GENOTYPE + DXDEP + HeadSize + nav + (years | RID)"))

fit_base <- brm(
  formula = formula_base,
  data = data,
  family = student(),
  prior = c(set_prior("normal(0, 10)", class = "b")),
  chains = 2,
  iter = 2000,
  warmup = 1000,
  refresh = 0,
  silent = 2
)

cat("--- Inizio ciclo sui volumi ---\n")
for (i in 1:length(volumes)) {
  name <- volumes[i]
  cat("Processing:", i, "/", length(volumes), "-", name, "\n")
  
  new_formula <- bf(as.formula(paste(name, "~ years*GENOTYPE + DXDEP + HeadSize + nav + (years | RID)")))
  
  tryCatch({
    model_bayesian <- update(fit_base, formula. = new_formula, newdata = data, refresh = 0)
    
    pred_frame <- info
    pred_frame$DXDEP <- factor(0, levels = levels(data$DXDEP))
    
    pred_frame$years <- 0
    y0 <- fitted(model_bayesian, newdata = pred_frame)[, "Estimate"]
    
    pred_frame$years <- 1
    y1 <- fitted(model_bayesian, newdata = pred_frame)[, "Estimate"]
    
    apc_values <- ((y1 - y0) / y0) * 100
    
    idx <- match(pred_frame$RID, newdata$RID)
    newdata[idx, name] <- apc_values
    
  }, error = function(e) { 
    cat("Error in volume:", name, "-", e$message, "\n") 
  })
}
newdata = newdata[, !(names(newdata) %in% c("nav", "HeadSize"))]
newdata4 = newdata






write.csv(newdata0, "newdata0.csv", row.names = FALSE)
write.csv(newdata1, "newdata1.csv", row.names = FALSE)
write.csv(newdata2, "newdata2.csv", row.names = FALSE)
write.csv(newdata3, "newdata3.csv", row.names = FALSE)
write.csv(newdata4, "newdata4.csv", row.names = FALSE)



# 5) LMM quantile regression <-- revealed some convergence issues -----------

datatot = read.csv("data.csv")
idno = which(datatot$nav != 1)
data = datatot[idno,]
data = data %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) == 1) %>%
  ungroup()

data$month = data$month / 12
colnames(data)[5] = "years"

volumes = colnames(data)[17:133]
u_rids = unique(data$RID)

info = data %>%
  group_by(RID) %>%
  summarise(
    DIAGNOSIS = first(DIAGNOSIS),
    GENOTYPE = first(GENOTYPE)
  )

newdata = data.frame(RID = u_rids)
newdata = left_join(newdata, info, by = "RID")
newdata[volumes] <- NA

for (name in volumes) {
  
  f = as.formula(paste(name, "~ years"))
  
  tryCatch({
    lqmm_mod = lqmm(fixed = f, 
                    random = ~ years, 
                    group = RID, 
                    tau = 0.5, 
                    data = data,
                    control = lqmmControl(method = "hashin"))
    
    f_eff = fixed.effects(lqmm_mod)
    re = ranef(lqmm_mod)
    
    slopes = f_eff["years"] + re[, "years"]
    intercepts = f_eff["(Intercept)"] + re[, "(Intercept)"]
    apc = (slopes / intercepts) * 100
    
    rid_names <- rownames(re)
    idx <- match(rid_names, newdata$RID)
    ok <- !is.na(idx)
    
    newdata[idx[ok], name] <- apc[ok]
  }, error = function(e) { cat("Errore nel volume:", name, "\n") })
}
























# Best model to estimate atrophy rates: compare with Wilcox-Test ------
datatot = read.csv("data.csv")
idno = which(datatot$nav != 1)
data = datatot[idno,]
data = data %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) == 1) %>%
  ungroup()

data$month = data$month / 12
colnames(data)[5] = "years"

data$DIAGNOSIS <- as.factor(data$DIAGNOSIS)
data$GENOTYPE <- as.factor(data$GENOTYPE)
data$DXDEP <- as.factor(data$DXDEP)

u_rids = unique(data$RID)
col_supra <- "Supratentorial.Total.Volume..Aseg...ST128SV."
hs_vector <- numeric(length(u_rids))
for (i in 1:length(u_rids)) {
  sub_p <- data[data$RID == u_rids[i], ]
  sub_p <- sub_p[order(sub_p$years), ]
  valid_vals <- sub_p[[col_supra]][!is.na(sub_p[[col_supra]])]
  hs_vector[i] <- if(length(valid_vals) > 0) valid_vals[1] else NA
}
hs_df <- data.frame(RID = u_rids, HeadSize = hs_vector)
data <- merge(data, hs_df, by = "RID", all.x = TRUE)

data_scaled <- data
data_scaled$HeadSize = data_scaled$HeadSize/1000 #transform to cm^3
target_vol <- colnames(data)[40] # CHOOSE HERE THE SIGNIFICANT VOLUME TO COMPARE

f1_scaled <- as.formula(paste(target_vol, "~ years + (years || RID)"))
f2_scaled <- as.formula(paste(target_vol, "~ years * GENOTYPE + (years | RID)"))
f3_scaled <- as.formula(paste(target_vol, "~ years * GENOTYPE + DXDEP + HeadSize + nav + (years | RID)"))

set.seed(123)
unique_ids <- unique(data_scaled$RID)
folds <- createFolds(unique_ids, k = 10)
mae_modA <- numeric(10)
mae_modB <- numeric(10)

for(i in 1:10) {
  test_ids <- unique_ids[folds[[i]]]
  train_data <- data_scaled[!(data_scaled$RID %in% test_ids), ]
  test_data <- data_scaled[data_scaled$RID %in% test_ids, ]
  
  tryCatch({
    mA_fold <- lmer(f2_scaled, data = train_data, REML = FALSE)
    mB_fold <- lmer(f3_scaled, data = train_data, REML = FALSE)
    
    predA <- predict(mA_fold, newdata = test_data, allow.new.levels = TRUE)
    predB <- predict(mB_fold, newdata = test_data, allow.new.levels = TRUE)
    
    mae_modA[i] <- mean(abs(predA - test_data[[target_vol]]), na.rm = TRUE)
    mae_modB[i] <- mean(abs(predB - test_data[[target_vol]]), na.rm = TRUE)
  }, error = function(e) {
    mae_modA[i] <- NA
    mae_modB[i] <- NA
  })
}

cv_comparison <- wilcox.test(mae_modA, mae_modB, paired = TRUE)

print(paste("Volume analizzato:", target_vol))
print(paste("MAE medio Modello A:", mean(mae_modA, na.rm = TRUE)))
print(paste("MAE medio Modello B:", mean(mae_modB, na.rm = TRUE)))
print(cv_comparison)

# Plot comparison --------------
uniforma_dati <- function(df) {
  df %>% mutate(
    DIAGNOSIS = as.factor(DIAGNOSIS),
    GENOTYPE = as.factor(GENOTYPE),
    RID = as.character(RID)
  )
}

newdata0 <- uniforma_dati(newdata0)
newdata1 <- uniforma_dati(newdata1)
newdata2 <- uniforma_dati(newdata2)
newdata3 <- uniforma_dati(newdata3)
newdata4 <- uniforma_dati(newdata4) 

cleandata0 <- uniforma_dati(cleandata0)
cleandata1 <- uniforma_dati(cleandata1)
cleandata2 <- uniforma_dati(cleandata2)
cleandata3 <- uniforma_dati(cleandata3)
cleandata4 <- uniforma_dati(cleandata4) 

newdata0$model <- "Quantile (Median)"
newdata1$model <- "LME Simple"
newdata2$model <- "LME + Covariates"
newdata3$model <- "LME + Covariates + etc"
newdata4$model <- "Bayesian LME"

cleandata0$model <- "Quantile (Median)"
cleandata1$model <- "LME Simple"
cleandata2$model <- "LME + Covariates"
cleandata3$model <- "LME + Covariates + etc"
cleandata4$model <- "Bayesian LME"

all_data <- bind_rows(cleandata0, cleandata1, cleandata4)

plot_data <- all_data %>%
  pivot_longer(cols = all_of(volumes), 
               names_to = "Volume", 
               values_to = "APC")

target_volumes <- volumes[1:min(1, length(volumes))]
target_volumes <- volumes[40]

df_filtered <- plot_data %>%
  filter(Volume %in% target_volumes)

ggplot(df_filtered, aes(x = APC, fill = model)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  facet_grid(model ~ DIAGNOSIS, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = paste("Absolute Counts of APC -", target_volumes),
       x = "Annual Percentage Change (%)",
       y = "Number of Patients")

ggplot(df_filtered, aes(x = DIAGNOSIS, y = APC, color = model)) +
  geom_jitter(alpha = 0.4, width = 0.2, size = 1) +
  geom_boxplot(alpha = 0.1, outlier.shape = NA, color = "black") +
  facet_wrap(~ model, nrow = 1) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Individual Patient Distribution per Model",
       x = "Diagnosis",
       y = "APC (%)")





# Outlier analysis -----------
## newdata0 --------
newdata0 = read.csv("newdata0.csv")
newdata_clean <- newdata0
all_cols <- colnames(newdata0)
volumes_cols <- all_cols[!(all_cols %in% c("RID", "DIAGNOSIS", "GENOTYPE"))]

ventricle_vars <- grep("Ventricle|CSF|Choroid|Plexus|Vessel|VentralDC|Inf.Lat.Vent", volumes_cols, value = TRUE, ignore.case = TRUE)
total_vol_vars <- grep("Total|Supratentorial|Cerebral.White.Matter|Gray.Matter.Volume|Subcortical.Volume", volumes_cols, value = TRUE, ignore.case = TRUE)
tissue_vars <- setdiff(volumes_cols, c(ventricle_vars, total_vol_vars))

for (v in tissue_vars) {
  valori <- suppressWarnings(as.numeric(as.character(newdata_clean[[v]])))
  newdata_clean[[v]][!is.na(valori) & (valori < -20 | valori > 5)] <- NA 
}

for (v in ventricle_vars) {
  valori <- suppressWarnings(as.numeric(as.character(newdata_clean[[v]])))
  newdata_clean[[v]][!is.na(valori) & (valori < -20 | valori > 100)] <- NA
}

for (v in total_vol_vars) {
  valori <- suppressWarnings(as.numeric(as.character(newdata_clean[[v]])))
  newdata_clean[[v]][!is.na(valori) & (valori < -5 | valori > 5)] <- NA
}
cleandata0 = newdata_clean
write.csv(cleandata0, "cleandata0.csv", row.names = FALSE)

na_counts <- colSums(is.na(cleandata0))
na_summary <- data.frame(
  Volume = names(na_counts),
  NA_Count = na_counts
) %>%
  filter(!(Volume %in% c("RID", "DIAGNOSIS", "GENOTYPE"))) %>%
  arrange(desc(NA_Count))
print(head(na_summary, 20))

## newdata1 -------------
newdata1 = read.csv("newdata1.csv")
newdata_clean <- newdata1
all_cols <- colnames(newdata1)
volumes_cols <- all_cols[!(all_cols %in% c("RID", "DIAGNOSIS", "GENOTYPE"))]

ventricle_vars <- grep("Ventricle|CSF|Choroid|Plexus|Vessel|VentralDC|Inf.Lat.Vent", volumes_cols, value = TRUE, ignore.case = TRUE)
total_vol_vars <- grep("Total|Supratentorial|Cerebral.White.Matter|Gray.Matter.Volume|Subcortical.Volume", volumes_cols, value = TRUE, ignore.case = TRUE)
tissue_vars <- setdiff(volumes_cols, c(ventricle_vars, total_vol_vars))

for (v in tissue_vars) {
  valori <- suppressWarnings(as.numeric(as.character(newdata_clean[[v]])))
  newdata_clean[[v]][!is.na(valori) & (valori < -20 | valori > 5)] <- NA
}

for (v in ventricle_vars) {
  valori <- suppressWarnings(as.numeric(as.character(newdata_clean[[v]])))
  newdata_clean[[v]][!is.na(valori) & (valori < -20 | valori > 100)] <- NA
}

for (v in total_vol_vars) {
  valori <- suppressWarnings(as.numeric(as.character(newdata_clean[[v]])))
  newdata_clean[[v]][!is.na(valori) & (valori < -5 | valori > 5)] <- NA
}
cleandata1 = newdata_clean


na_counts <- colSums(is.na(cleandata1))
na_summary <- data.frame(
  Volume = names(na_counts),
  NA_Count = na_counts
) %>%
  filter(!(Volume %in% c("RID", "DIAGNOSIS", "GENOTYPE"))) %>%
  arrange(desc(NA_Count))
print(head(na_summary, 20))

cleandata1 = cleandata1[, !(names(newdata) %in% c("Supratentorial.Total.Volume..Aseg...ST128SV."))]
write.csv(cleandata1, "cleandata1.csv", row.names = FALSE)

## newdata2 -----------
newdata2 = read.csv("newdata2.csv")
newdata_clean <- newdata2
all_cols <- colnames(newdata2)
volumes_cols <- all_cols[!(all_cols %in% c("RID", "DIAGNOSIS", "GENOTYPE"))]

ventricle_vars <- grep("Ventricle|CSF|Choroid|Plexus|Vessel|VentralDC|Inf.Lat.Vent", volumes_cols, value = TRUE, ignore.case = TRUE)
total_vol_vars <- grep("Total|Supratentorial|Cerebral.White.Matter|Gray.Matter.Volume|Subcortical.Volume", volumes_cols, value = TRUE, ignore.case = TRUE)
tissue_vars <- setdiff(volumes_cols, c(ventricle_vars, total_vol_vars))

for (v in tissue_vars) {
  valori <- suppressWarnings(as.numeric(as.character(newdata_clean[[v]])))
  newdata_clean[[v]][!is.na(valori) & (valori < -20 | valori > 5)] <- NA
}

for (v in ventricle_vars) {
  valori <- suppressWarnings(as.numeric(as.character(newdata_clean[[v]])))
  newdata_clean[[v]][!is.na(valori) & (valori < -20 | valori > 100)] <- NA
}

for (v in total_vol_vars) {
  valori <- suppressWarnings(as.numeric(as.character(newdata_clean[[v]])))
  newdata_clean[[v]][!is.na(valori) & (valori < -5 | valori > 5)] <- NA
}
cleandata2 = newdata_clean

na_counts <- colSums(is.na(cleandata2))
na_summary <- data.frame(
  Volume = names(na_counts),
  NA_Count = na_counts
) %>%
  filter(!(Volume %in% c("RID", "DIAGNOSIS", "GENOTYPE"))) %>%
  arrange(desc(NA_Count))
print(head(na_summary, 20))

cleandata2 = cleandata2[, !(names(newdata) %in% c("Supratentorial.Total.Volume..Aseg...ST128SV."))]
write.csv(cleandata2, "cleandata2.csv", row.names = FALSE)

## newdata3 ----------
newdata3 = read.csv("newdata3.csv")
newdata_clean <- newdata3
all_cols <- colnames(newdata3)
volumes_cols <- all_cols[!(all_cols %in% c("RID", "DIAGNOSIS", "GENOTYPE"))]

ventricle_vars <- grep("Ventricle|CSF|Choroid|Plexus|Vessel|VentralDC|Inf.Lat.Vent", volumes_cols, value = TRUE, ignore.case = TRUE)
total_vol_vars <- grep("Total|Supratentorial|Cerebral.White.Matter|Gray.Matter.Volume|Subcortical.Volume", volumes_cols, value = TRUE, ignore.case = TRUE)
tissue_vars <- setdiff(volumes_cols, c(ventricle_vars, total_vol_vars))

for (v in tissue_vars) {
  valori <- suppressWarnings(as.numeric(as.character(newdata_clean[[v]])))
  newdata_clean[[v]][!is.na(valori) & (valori < -20 | valori > 5)] <- NA
}

for (v in ventricle_vars) {
  valori <- suppressWarnings(as.numeric(as.character(newdata_clean[[v]])))
  newdata_clean[[v]][!is.na(valori) & (valori < -20 | valori > 100)] <- NA
}

for (v in total_vol_vars) {
  valori <- suppressWarnings(as.numeric(as.character(newdata_clean[[v]])))
  newdata_clean[[v]][!is.na(valori) & (valori < -5 | valori > 5)] <- NA
}
cleandata3 = newdata_clean

na_counts <- colSums(is.na(cleandata3))
na_summary <- data.frame(
  Volume = names(na_counts),
  NA_Count = na_counts
) %>%
  filter(!(Volume %in% c("RID", "DIAGNOSIS", "GENOTYPE"))) %>%
  arrange(desc(NA_Count))
print(head(na_summary, 20))

cleandata3 = cleandata3[, !(names(newdata) %in% c("Supratentorial.Total.Volume..Aseg...ST128SV."))]

write.csv(cleandata3, "cleandata3.csv", row.names = FALSE)
## newdata4 -------------
newdata4 = read.csv("newdata4.csv")
newdata_clean <- newdata4
all_cols <- colnames(newdata4)
volumes_cols <- all_cols[!(all_cols %in% c("RID", "DIAGNOSIS", "GENOTYPE"))]

ventricle_vars <- grep("Ventricle|CSF|Choroid|Plexus|Vessel|VentralDC|Inf.Lat.Vent", volumes_cols, value = TRUE, ignore.case = TRUE)
total_vol_vars <- grep("Total|Supratentorial|Cerebral.White.Matter|Gray.Matter.Volume|Subcortical.Volume", volumes_cols, value = TRUE, ignore.case = TRUE)
tissue_vars <- setdiff(volumes_cols, c(ventricle_vars, total_vol_vars))

for (v in tissue_vars) {
  valori <- suppressWarnings(as.numeric(as.character(newdata_clean[[v]])))
  newdata_clean[[v]][!is.na(valori) & (valori < -20 | valori > 5)] <- NA
}

for (v in ventricle_vars) {
  valori <- suppressWarnings(as.numeric(as.character(newdata_clean[[v]])))
  newdata_clean[[v]][!is.na(valori) & (valori < -20 | valori > 100)] <- NA
}

for (v in total_vol_vars) {
  valori <- suppressWarnings(as.numeric(as.character(newdata_clean[[v]])))
  newdata_clean[[v]][!is.na(valori) & (valori < -5 | valori > 5)] <- NA
}
cleandata4 = newdata_clean
write.csv(cleandata4, "cleandata4.csv", row.names = FALSE)

na_counts <- colSums(is.na(cleandata4))
na_summary <- data.frame(
  Volume = names(na_counts),
  NA_Count = na_counts
) %>%
  filter(!(Volume %in% c("RID", "DIAGNOSIS", "GENOTYPE"))) %>%
  arrange(desc(NA_Count))
print(head(na_summary, 20))




