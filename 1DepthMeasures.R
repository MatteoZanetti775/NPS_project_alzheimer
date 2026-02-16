library(MASS) 
library(dplyr)
library(rgl) 
library(DepthProc) 
library(hexbin) 
library(aplpack) 
library(robustbase)
library(fda.usc)
library(fda)
set.seed(20)


data = read.csv("data.csv")
dim(data)
colnames(data)


id = which(data$nav==1)
sub = data[id, 17:25] #Let's take only 1-visit patients only
dim(sub)
id1 = which(data[id, ]$DIAGNOSIS ==1)
id2 = which(data[id, ]$DIAGNOSIS ==2)
id4 = which(data[id, ]$DIAGNOSIS ==4)
sub1 = sub[id1,]
sub2 = sub[id2,]
sub4 = sub[id4,]

#Depth measures ----------
tukey_depth1 = depth(u = sub1, method='Tukey')
tukey_depth2 = depth(u = sub2, method='Tukey')
tukey_depth4 = depth(u = sub4, method='Tukey')
sub[which.max(tukey_depth1),]
sub[which.max(tukey_depth2),]
sub[which.max(tukey_depth4),]
#s_depth = mrfDepth::sdepth(data)$depthZ #This depth works only for p<=2
proj_depth1 = depth(u = sub1, method='Projection')
proj_depth2 = depth(u = sub2, method='Projection')
proj_depth4 = depth(u = sub4, method='Projection')
maha_depth1 = depth(sub1, method ='Mahalanobis')
maha_depth2 = depth(sub2, method ='Mahalanobis')
maha_depth4 = depth(sub4, method ='Mahalanobis')
#Pairwise bagplots => not a good idea to identify outliers
bagplot1 = with(sub1, aplpack::bagplot(Right.Pars.Opercularis.Volume..ST104CV., Right.Posterior.Cingulate.Volume..ST109CV.)) 
bagplot2 = with(sub2, aplpack::bagplot(Right.Pars.Opercularis.Volume..ST104CV., Right.Posterior.Cingulate.Volume..ST109CV.)) 
bagplot4 = with(sub4, aplpack::bagplot(Right.Pars.Opercularis.Volume..ST104CV., Right.Posterior.Cingulate.Volume..ST109CV.)) 
bagplotpairs = with(sub, aplpack::bagplot.pairs(sub1)) 
bagplotpairs = with(sub, aplpack::bagplot.pairs(sub2)) 
bagplotpairs = with(sub, aplpack::bagplot.pairs(sub4)) 
out = bagplot1$pxy.outlier #doesn't work for bagplotpairs
id_out_obs = which(apply(sub1, 1, function(x) all(x %in% out)))
id_out_obs
id_out_obs = which(apply(sub2, 1, function(x) all(x %in% out)))
id_out_obs
id_out_obs = which(apply(sub4, 1, function(x) all(x %in% out)))
id_out_obs
#Plot(ID, depth)
plot(data[id1, "RID"], tukey_depth1) #Majority of Turkey depths are compressed to 0 (they are not necessarily outliers)
plot(data[id1, "RID"], tukey_depth1, log="y")
plot(rank(tukey_depth1), tukey_depth1)
plot(data[id2, "RID"], tukey_depth2, log="y")
plot(rank(tukey_depth2), tukey_depth2)
plot(data[id4, "RID"], tukey_depth4, log="y")
plot(rank(tukey_depth4), tukey_depth4)
plot(data[id1, "RID"], maha_depth1)
plot(data[id2, "RID"], maha_depth2)
plot(data[id4, "RID"], maha_depth4)
plot(rank(maha_depth1), maha_depth1) #Cannot really identify any elbow
plot(rank(maha_depth2), maha_depth2)
plot(rank(maha_depth4), maha_depth4)
plot(data[id1, "RID"], proj_depth1)
plot(data[id2, "RID"], proj_depth2)
plot(data[id4, "RID"], proj_depth4)
plot(rank(proj_depth1), proj_depth1) #Cannot really identify any elbow
plot(rank(proj_depth2), proj_depth2) 
plot(rank(proj_depth4), proj_depth4) 
#Hist(depth): depth concentered towards 0 due to d>>, highly skewed dsitr
hist(tukey_depth1)
hist(tukey_depth2)
hist(tukey_depth4)
hist(maha_depth1) #There seem to be a threshold for outliers at 0.05
hist(maha_depth2)
hist(maha_depth4)
#QQplot(tureky,maha): Tukey is collapsing --> a lot of values ≈ 0 => not really able to order pts  
qqplot(sort(tukey_depth1), sort(maha_depth1))
qqplot(sort(tukey_depth2), sort(maha_depth2))
qqplot(sort(tukey_depth4), sort(maha_depth4))
#DDplot(turkey, maha): the 2 depths do not really agree
plot(tukey_depth1, maha_depth1) 
#Outlier detection using quantiles
threshold = quantile(proj_depth1, 0.05)
outliers = which(proj_depth1 <= threshold)
data[id1[outliers], "RID"]
threshold = quantile(proj_depth2, 0.05)
outliers = which(proj_depth2 <= threshold)
data[id2[outliers], "RID"]
threshold = quantile(proj_depth4, 0.05)
outliers = which(proj_depth4 <= threshold)
data[id4[outliers], "RID"]
#Outlier detection using median, median absolute deviation
median_depth = median(proj_depth1)
mad_depth = mad(proj_depth1)
outliers = which(proj_depth1 < median_depth - 3*mad_depth)
outliers
data[id1[outliers], "RID"]
median_depth = median(proj_depth2)
mad_depth = mad(proj_depth2)
outliers = which(proj_depth2 < median_depth - 3*mad_depth)
outliers
data[id2[outliers], "RID"]
median_depth = median(proj_depth4)
mad_depth = mad(proj_depth4)
outliers = which(proj_depth4 < median_depth - 3*mad_depth)
outliers
data[id4[outliers], "RID"]


#STANDARDIZATIONS of the VOLS ------------
datastd = data
for(i in 17:133){
  datastd[,i] <- (datastd[,i] - mean(datastd[,i], na.rm = TRUE)) / sd(datastd[,i], na.rm = TRUE)
}
#Creation of visit-datasets ----------
id9 = which(data$nav == 9)
id8 = which(data$nav == 8)
id7 = which(data$nav == 7)
id6 = which(data$nav == 6)
id5 = which(data$nav == 5)
id4 = which(data$nav == 4)
id3 = which(data$nav == 3)
id2 = which(data$nav == 2)
data9visits = data[id9,c(1,7,9,17:133)]
eight = data[id8,c(1,7,9,17:133)]
seven = data[id7,c(1,7,9,17:133)]
six = data[id6,c(1,7,9,17:133)]
five = data[id5,c(1,7,9,17:133)]
four = data[id4,c(1,7,9,17:133)]
three = data[id3,c(1,7,9,17:133)]
two = data[id2,c(1,7,9,17:133)]



stable <- ungroup(filter(group_by(data, RID), n_distinct(DIAGNOSIS) == 1))
stable1 <- stable[stable$DIAGNOSIS == 1, ]
stable2 <- stable[stable$DIAGNOSIS == 2, ]
stable4 <- stable[stable$DIAGNOSIS == 4, ]
#PLOT WITHOUT SMOOTHING ------------
ggplot(data, aes(x = age, y =  Right.Hippocampus.Volume..ST77SV., color = factor(DIAGNOSIS), group = RID)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(color = "Diagnosis", x = "Age", y = "Right Hippocampus Volume", title = "All patients")
ggplot(data9visits, aes(x = age, y = Right.Hippocampus.Volume..ST77SV., color = factor(DIAGNOSIS), group = RID)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(color = "Diagnosis", x = "Age", y = "Right Hippocampus Volume", title = "9-visits only patients")
ggplot(eight, aes(x = age, y = Right.Hippocampus.Volume..ST77SV., color = factor(DIAGNOSIS), group = RID)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(color = "Diagnosis", x = "Age", y = "Right Hippocampus Volume", title = "8-visits only patients")
ggplot(seven, aes(x = age, y = Right.Hippocampus.Volume..ST77SV., color = factor(DIAGNOSIS), group = RID)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(color = "Diagnosis", x = "Age", y = "Right Hippocampus Volume", title = "7-visits only patients")
ggplot(six, aes(x = age, y = Right.Hippocampus.Volume..ST77SV., color = factor(DIAGNOSIS), group = RID)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(color = "Diagnosis", x = "Age", y = "Right Hippocampus Volume", title = "6-visits only patients")
ggplot(five, aes(x = age, y = Right.Hippocampus.Volume..ST77SV., color = factor(DIAGNOSIS), group = RID)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(color = "Diagnosis", x = "Age", y = "Right Hippocampus Volume", title = "5-visits only patients")
ggplot(four, aes(x = age, y = Right.Hippocampus.Volume..ST77SV., color = factor(DIAGNOSIS), group = RID)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(color = "Diagnosis", x = "Age", y = "Right Hippocampus Volume", title = "4-visits only patients")
ggplot(three, aes(x = age, y = Right.Hippocampus.Volume..ST77SV., color = factor(DIAGNOSIS), group = RID)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(color = "Diagnosis", x = "Age", y = "Right Hippocampus Volume", title = "3-visits only patients")
ggplot(two, aes(x = age, y = Right.Hippocampus.Volume..ST77SV., color = factor(DIAGNOSIS), group = RID)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(color = "Diagnosis", x = "Age", y = "Right Hippocampus Volume", title = "2-visits only patients")
ggplot(stable1, aes(x = age, y = Right.Hippocampus.Volume..ST77SV., color = factor(DIAGNOSIS), group = RID)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(color = "Diagnosis", x = "Age", y = "Right Hippocampus Volume", title = "patients with stable diagnosis 1")
ggplot(stable2, aes(x = age, y = Right.Hippocampus.Volume..ST77SV., color = factor(DIAGNOSIS), group = RID)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(color = "Diagnosis", x = "Age", y = "Right Hippocampus Volume", title = "patients with stable diagnosis 2")
ggplot(stable4, aes(x = age, y = Right.Hippocampus.Volume..ST77SV., color = factor(DIAGNOSIS), group = RID)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(color = "Diagnosis", x = "Age", y = "Right Hippocampus Volume", title = "patients with stable diagnosis 4")
#PLOT WITH SMOOTHING -------
library(fda)
library(dplyr)

smooth_volume_plot <- function(df, volume_col, plot_title=NULL) {
  curve_list <- df %>% 
    group_by(RID) %>% 
    summarise(
      age = list(age),
      volume = list(!!sym(volume_col))
    )
  
  fd_list <- list()
  age_grid <- seq(min(df$age), max(df$age), length.out = 200)
  
  # Plot vuoto
  plot(NULL, xlim=range(age_grid), ylim=range(df[[volume_col]], na.rm=TRUE),
       xlab="Age", ylab=volume_col, main=plot_title)
  grid(nx=NA, ny=NULL, col="lightgray", lty="dotted")
  
  colors <- rainbow(nrow(curve_list))
  
  for(i in 1:nrow(curve_list)){
    ages <- curve_list$age[[i]]
    vols <- curve_list$volume[[i]]
    n_obs <- length(ages)
    
    
    if(any(vols <= 0, na.rm = TRUE)){
      lines(ages, vols, col=colors[i], lwd=2)
      next
    }
    
    
    if(n_obs <= 5){
      lines(ages, vols, col=colors[i], lwd=2)
      next
    }
    
    
    vols_log <- log(vols)
    nbasis <- min(n_obs, 6)
    norder <- ifelse(nbasis < 4, nbasis, 4)
    basis <- create.bspline.basis(rangeval = range(ages), nbasis = nbasis, norder = norder)
    
    Lfdobj <- ifelse(n_obs < 6, 1, 2)
    lambda <- ifelse(n_obs < 6, 1e-6, 1e-2)
    fdParobj <- fdPar(basis, Lfdobj = Lfdobj, lambda = lambda)
    
    fd <- smooth.basis(ages, vols_log, fdParobj)$fd
    fd_list[[i]] <- fd
    
    r <- fd$basis$rangeval
    idx <- which(age_grid >= r[1] & age_grid <= r[2])
    y <- exp(eval.fd(age_grid[idx], fd))  # Torna alla scala originale
    lines(age_grid[idx], y, col=colors[i], lwd=2)
  }
  
  return(fd_list)
}


fd9 <- smooth_volume_plot(data9visits, "Right.Hippocampus.Volume..ST77SV.", "9-visits")
fd8 <- smooth_volume_plot(eight, "Right.Hippocampus.Volume..ST77SV.", "8-visits")
fd7 <- smooth_volume_plot(seven, "Right.Hippocampus.Volume..ST77SV.", "7-visits")
fd6 <- smooth_volume_plot(six, "Right.Hippocampus.Volume..ST77SV.", "6-visits")
fd5 <- smooth_volume_plot(five, "Right.Hippocampus.Volume..ST77SV.", "5-visits")
fd4 <- smooth_volume_plot(four, "Right.Hippocampus.Volume..ST77SV.", "4-visits")
fd3 <- smooth_volume_plot(three, "Right.Hippocampus.Volume..ST77SV.", "3-visits")
fd2 <- smooth_volume_plot(two, "Right.Hippocampus.Volume..ST77SV.", "2-visits")




#FUNCTIONAL DEPTH MEASURES: Modified Band Depth multivariate --------
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
  
  
write.csv(depth_merged, "depth_stable_patients.csv", row.names = FALSE)
