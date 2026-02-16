#This file uses atrophy rates computed in step 3
#Can use cleandata with outliers too
library(dplyr)
library(coin)
library(boot)
library(tidyr)
library(effsize)
library(ggcorrplot)

datatot <- read.csv("data.csv")
newdata = read.csv("cleandata0.csv") 
sum(is.na(newdata))

data <- datatot %>%
  filter(nav != 1) %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) == 1) %>%
  ungroup()

data$month <- data$month / 12
colnames(data)[5] <- "years"
data$DIAGNOSIS <- as.factor(data$DIAGNOSIS)

volumes <- colnames(data)[17:133]
u_rids <- unique(data$RID)

results_perm <- data.frame(
  Volume = volumes,
  p_value = NA,
  stat = NA
)

# Univariate Permutational Tests ---------
for (i in 1:length(volumes)) {
  vol_name <- volumes[i]
  f <- as.formula(paste(vol_name, "~ DIAGNOSIS"))
  
  tryCatch({
    test_res <- independence_test(f, 
                                  data = newdata, 
                                  distribution = approximate(nresample = 10000))
    results_perm$p_value[i] <- pvalue(test_res)
    results_perm$stat[i] <- statistic(test_res)
  }, error = function(e) {})
}

results_perm$p_adj <- p.adjust(results_perm$p_value, method = "BH")
sig_volumes <- results_perm %>% filter(p_adj < 0.05) %>% pull(Volume)
non_sig_volumes <- volumes[!(volumes %in% sig_volumes)]
non_sig_volumes

# Effect anlaysis: interested not only in the p-value but also to quantify the difference ------
significant_vols_list <- results_perm %>% filter(p_adj < 0.05) %>% pull(Volume)

full_summary_table <- data.frame()

for (vol in significant_vols_list) {
  diag_levels <- unique(newdata$DIAGNOSIS)
  
  for (i in 1:length(diag_levels)) {
    d_label <- diag_levels[i]
    sub_data <- newdata[newdata$DIAGNOSIS == d_label, ]
    
    boot_fn <- function(d, indices) {
      return(mean(d[indices, vol], na.rm = TRUE))
    }
    
    b_obj <- boot::boot(data = sub_data, statistic = boot_fn, R = 1000)
    b_ci <- boot::boot.ci(b_obj, type = "perc")
    
    full_summary_table <- rbind(full_summary_table, data.frame(
      Volume = vol,
      Diagnosis = d_label,
      Mean_APC = b_obj$t0,
      CI_lower = b_ci$percent[4],
      CI_upper = b_ci$percent[5]
    ))
  }
}

effect_sizes_base <- full_summary_table %>%
  filter(Diagnosis %in% c("1", "4")) %>%
  pivot_wider(names_from = Diagnosis, values_from = c(Mean_APC, CI_lower, CI_upper)) %>%
  mutate(Atrophy_Delta = Mean_APC_4 - Mean_APC_1)

effect_sizes_atrophy <- effect_sizes_base %>%
  filter(Atrophy_Delta < 0) %>%
  arrange(Atrophy_Delta)

effect_sizes_expansion <- effect_sizes_base %>%
  filter(Atrophy_Delta > 0) %>%
  arrange(desc(Atrophy_Delta))

top_atrophy <- head(effect_sizes_atrophy, 3)
top_espansion <- head(effect_sizes_expansion, 2)

plot_combined <- rbind(top_atrophy, top_espansion)
top_brain_markers_final <- plot_combined$Volume

print(plot_combined %>% select(Volume, Mean_APC_1, Mean_APC_4, Atrophy_Delta))



ggplot(plot_combined, aes(x = reorder(Volume, Atrophy_Delta), y = Atrophy_Delta)) +
  geom_segment(aes(xend = Volume, yend = 0), color = "grey80") +
  geom_point(aes(color = Atrophy_Delta > 0), size = 4) +
  scale_color_manual(values = c("blue", "red"), labels = c("Atrophy", "Expansion")) +
  coord_flip() +
  theme_minimal() +
  labs(title = "AD-Specific Structural Changes Ranking",
       subtitle = "Top Atrophy and Expansion Markers (AD vs CN Delta)",
       x = "Brain Region",
       y = "Delta APC (%)",
       color = "Change Type")


# Some plots -------------
delta_boot_results <- data.frame()

for (vol in top_brain_markers_final) {
  
  sub_data <- newdata[newdata$DIAGNOSIS %in% c("1","4"), ]
  sub_data <- sub_data[!is.na(sub_data[[vol]]), ]
  
  boot_diff_fn <- function(d, indices) {
    d2 <- d[indices, ]
    mean(d2[d2$DIAGNOSIS=="4", vol], na.rm=TRUE) -
      mean(d2[d2$DIAGNOSIS=="1", vol], na.rm=TRUE)
  }
  
  b_obj <- boot(data = sub_data, statistic = boot_diff_fn, R = 1000)
  b_ci <- boot.ci(b_obj, type="perc")
  
  delta_boot_results <- rbind(delta_boot_results,
                              data.frame(
                                Volume = vol,
                                Delta = b_obj$t0,
                                CI_lower = b_ci$percent[4],
                                CI_upper = b_ci$percent[5]
                              ))
}

volume_labels <- c(
  "Left.Frontal.Pole.Volume..ST25CV." = "Left Frontal Pole",
  "Left.Thalamus.Volume..ST11SV." = "Left Thalamus",
  "Left.Entorhinal.Volume..ST24CV." = "Left Entorhinal Cortex",
  "Right.Vessel.Volume..ST125SV." = "Right Vessel",
  "Left.Inferior.Lateral.Ventricle.Volume..ST30SV." = "Left Inferior Lateral Ventricle"
)


trend_data_plot <- full_summary_table %>%
  filter(Volume %in% top_brain_markers_final) %>%
  mutate(
    ClinicalGroup = factor(Diagnosis, levels = c("1","2","4"), labels = c("CN","MCI","AD")),
    VolumeLabel = volume_labels[Volume]
  )

plot_means <- ggplot(trend_data_plot,
                     aes(x = ClinicalGroup,
                         y = Mean_APC,
                         group = VolumeLabel,
                         color = VolumeLabel)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.15, linewidth = 0.8) +
  scale_color_manual(values = c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e")) +
  theme_minimal(base_size = 16) +
  labs(title = "Bootstrap Mean Annual Atrophy Rate",
       x = "Clinical Group",
       y = "Mean APC (%)",
       color = "Brain Region")

plot_means

plot_delta <- ggplot(delta_boot_results,
                     aes(x = reorder(Volume, Delta),
                         y = Delta)) +
  geom_point(aes(color = Delta > 0), size = 4) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.2, linewidth = 0.8) +
  scale_color_manual(values = c("blue","red"), labels = c("Atrophy","Expansion")) +
  scale_x_discrete(labels = volume_labels) +
  coord_flip() +
  theme_minimal(base_size = 16) +
  labs(title = "Bootstrap CI of AD–CN Differences",
       x = "Brain Region",
       y = "Delta APC (%)",
       color = "Change Type")

plot_delta







# Top 15 significant volumes --------
indices_in_volumes_vector <- match(top_brain_markers_final, volumes)
print(indices_in_volumes_vector)
print(volumes[indices_in_volumes_vector])

colSums(is.na(newdata[, (names(newdata) %in%  volumes[indices_in_volumes_vector])]))









# Overlap between groups --------
separation_results <- data.frame()

for (vol in top_brain_markers_final) {
  sub_data <- newdata[newdata$DIAGNOSIS %in% c("1", "4"), ]
  sub_data <- sub_data[!is.na(sub_data[[vol]]), ]
  
  cohen <- cohen.d(as.formula(paste(vol, "~ DIAGNOSIS")), data = sub_data)
  
  separation_results <- rbind(separation_results, data.frame(
    Volume = vol,
    Cohen_d = cohen$estimate,
    Magnitude = cohen$magnitude
  ))
}

separation_results <- separation_results[order(abs(separation_results$Cohen_d), decreasing = TRUE), ]
print(separation_results)


# Progression stage ---------
trend_data <- full_summary_table %>%
  filter(Volume %in% head(top_atrophy$Volume, 5))

trend_data_plot <- trend_data %>%
  mutate(ClinicalGroup = factor(Diagnosis, levels = c("1","2","4"),
                                labels = c("CN","MCI","AD")))

ggplot(trend_data_plot, aes(x = ClinicalGroup, y = Mean_APC, group = Volume, color = Volume)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.1) +
  scale_color_manual(values = c(
    "Left.Frontal.Pole.Volume..ST25CV." = "#1b9e77",
    "Left.Thalamus.Volume..ST11SV." = "#d95f02",
    "Left.Entorhinal.Volume..ST24CV." = "#7570b3",
    "Right.Vessel.Volume..ST125SV." = "#e7298a",
    "Left.Inferior.Lateral.Ventricle.Volume..ST30SV." = "#66a61e"
  ),
  labels = c(
    "Left.Frontal.Pole.Volume..ST25CV." = "Left Frontal Pole",
    "Left.Thalamus.Volume..ST11SV." = "Left Thalamus",
    "Left.Entorhinal.Volume..ST24CV." = "Left Entorhinal Cortex",
    "Right.Vessel.Volume..ST125SV." = "Right Vessel",
    "Left.Inferior.Lateral.Ventricle.Volume..ST30SV." = "Left Inferior Lateral Ventricle"
  )) +
  theme_minimal() +
  labs(title = "Progression Profile: Atrophy across Disease Stages",
       x = "Clinical Group",
       y = "Mean Annual Percentage Change (%)",
       color = "Brain Region")

#confidence intervals do not overlap => statistically significant difference


# Correlation between top markers
top_data <- newdata[, top_brain_markers_final]
corr_matrix <- cor(top_data, use = "pairwise.complete.obs")

ggcorrplot(corr_matrix, 
           hc.order = TRUE, 
           type = "lower", 
           lab = TRUE, 
           title = "Redundancy Analysis: Correlation between Top Markers")


# GENOTYPE IMPACT -----------
newdata$GENOTYPE <- as.factor(newdata$GENOTYPE)
genotype_impact <- data.frame()
for (vol in head(top_atrophy$Volume, 5)) {
  f_gen <- as.formula(paste(vol, "~ GENOTYPE"))
  test_gen <- independence_test(f_gen, data = newdata[newdata$DIAGNOSIS == "4", ],
                                distribution = approximate(nresample = 5000))
  
  genotype_impact <- rbind(genotype_impact, data.frame(
    Volume = vol,
    Genotype_p_val = pvalue(test_gen)
  ))
}
print(genotype_impact)



  