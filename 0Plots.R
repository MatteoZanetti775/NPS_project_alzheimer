library(dplyr)
library(data.table) 
library(ggplot2)
data = read.csv("data.csv")
#Add a col which takes into account the number of visits AVAILABLE for each patient
data$nav <- ave(data$RID, data$RID, FUN = length)
data <- relocate(data, nav, .after = entry)
dim(data)
colnames(data)
length(unique(data$RID))
table(data$GENOTYPE)




# 0. distribution of research group overall - not useful -------
distribuzione_categoria <- data %>%
  group_by(RID, entry_research_group) %>%
  summarise(.groups = "drop_last") %>%
  summarise(n_categorie = n(), .groups = "drop")
table(data$entry_research_group)
ggplot(data, aes(x = entry_research_group)) +
  geom_bar(fill = "steelblue") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  labs(
    title = "Distribution of diagnosed research group",
    x = "Levels",
    y = "Counts"
  ) +
  theme_minimal(base_size = 14)
## 0.2 Distribution of entry research group per patient -------
categoria_per_paziente <- data %>%
  distinct(RID, entry_research_group)
table(categoria_per_paziente$entry_research_group)
ggplot(categoria_per_paziente, aes(x = entry_research_group)) +
  geom_bar(fill = "steelblue") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5)+
  labs(
    title = "Distribuzione dei pazienti per entry research group",
    x = "Entry research group",
    y = "Numero di pazienti"
  ) +
  theme_minimal(base_size = 14)


# 0.3 Time-lines different DIAGNOSIS each patient has ----------
patients_change <- data %>%
  group_by(RID) %>%
  summarise(
    diag_min = min(DIAGNOSIS),
    diag_max = max(DIAGNOSIS),
    diag_seq = paste(unique(DIAGNOSIS), collapse = "-"),
    .groups = "drop"
  )
patients_regress <- patients_change %>%
  filter(diag_min < diag_max & diag_seq %in% c("2-1","3-1","3-2-1"))
patients_1_3 <- patients_change %>%
  filter(diag_seq == "1-3") %>% pull(RID)
patients_2_3 <- patients_change %>%
  filter(diag_seq == "2-3") %>% pull(RID)
patients_1_2_3 <- patients_change %>%
  filter(diag_seq == "1-2-3") %>% pull(RID)

data_1_2_3 <- data %>% filter(RID %in% patients_1_2_3)
data_1_3 <- data %>% filter(RID %in% patients_1_3)
data_2_3 <- data %>% filter(RID %in% patients_2_3)

# Plot 1 -> 2 -> 3
ggplot(data_1_2_3, aes(x = month, y = factor(RID), group = RID, color = factor(DIAGNOSIS))) +
  geom_line(alpha = 0.5) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "Patients with diagnosis change 1 -> 2 -> 3",
    x = "Month from baseline",
    y = "Patient ID",
    color = "Diagnosis"
  ) +
  theme_minimal(base_size = 12)
# Plot 1 -> 3
ggplot(data_1_3, aes(x = month, y = factor(RID), group = RID, color = factor(DIAGNOSIS))) +
  geom_line(alpha = 0.5) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "Patients with diagnosis change 1 -> 3",
    x = "Month from baseline",
    y = "Patient ID",
    color = "Diagnosis"
  ) +
  theme_minimal(base_size = 12)
# Plot 2 -> 3
ggplot(data_2_3, aes(x = month, y = factor(RID), group = RID, color = factor(DIAGNOSIS))) +
  geom_line(alpha = 0.5) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "Patients with diagnosis change 2 -> 3",
    x = "Month from baseline",
    y = "Patient ID",
    color = "Diagnosis"
  ) +
  theme_minimal(base_size = 12)

patients_regress # non ci sono pazienti che migliorano nel tempo

## Istograms for changed DIAGNOSIS --------
patients_change <- patients_change %>%
  mutate(
    progression_category = case_when(
      diag_seq == "1-3"     ~ "1 -> 3",
      diag_seq == "2-3"     ~ "2 -> 3",
      diag_seq == "1-2-3"   ~ "1 -> 2 -> 3",
      TRUE                   ~ NA_character_
    )
  )
patients_progression <- patients_change %>%
  filter(!is.na(progression_category)) %>%
  left_join(data %>% distinct(RID, entry_research_group), by = "RID")
progression_counts <- patients_progression %>%
  group_by(progression_category, entry_research_group) %>%
  summarise(n_patients = n(), .groups = "drop")

ggplot(progression_counts, aes(x = progression_category, y = n_patients, fill = entry_research_group)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = n_patients), 
            position = position_stack(vjust = 0.5), size = 4) +
  labs(
    title = "Number of patients by diagnosis progression and baseline group",
    x = "Progression category",
    y = "Number of patients",
    fill = "Entry research group"
  ) +
  theme_minimal(base_size = 14)

### Dig deeper for those 3 individuals that should start from 2 but they don't: INTERESTING CASES --------
groups <- c("AD", "CN", "SMC")
data_2_3_subset <- data %>%
  filter(RID %in% patients_2_3, entry_research_group %in% groups)
ggplot(data_2_3_subset, aes(x = month, y = factor(RID), group = RID, color = factor(DIAGNOSIS))) +
  geom_line(alpha = 0.5) +
  geom_point(size = 2, alpha = 0.7) +
  facet_wrap(~ entry_research_group, ncol = 1, scales = "free_y") +  # asse x comune, y libera
  labs(
    title = "Timeline of visits for patients progressing 2 -> 3",
    x = "Month from baseline",
    y = "Patient ID",
    color = "Diagnosis"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Istograms for NOT changed DIAGNOSIS --------
patients_stable <- patients_change %>%
  filter(diag_min == diag_max) %>%
  mutate(
    stable_category = case_when(
      diag_min == 1 ~ "1 -> 1",
      diag_min == 2 ~ "2 -> 2",
      diag_min == 3 ~ "3 -> 3",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(stable_category)) %>%
  left_join(data %>% distinct(RID, entry_research_group), by = "RID")
stable_counts <- patients_stable %>%
  group_by(stable_category, entry_research_group) %>%
  summarise(n_patients = n(), .groups = "drop")
ggplot(stable_counts, aes(x = stable_category, y = n_patients, fill = entry_research_group)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = n_patients), position = position_stack(vjust = 0.5), size = 4) +
  labs(
    title = "Number of patients with stable diagnosis by baseline group",
    x = "Stable diagnosis category",
    y = "Number of patients",
    fill = "Entry research group"
  ) +
  theme_minimal(base_size = 14)
### Dig deeper for suspicious cases ---------
stable_pivot <- patients_stable %>%
  group_by(stable_category, entry_research_group) %>%
  summarise(n_patients = n(), .groups = "drop") %>%
  pivot_wider(names_from = entry_research_group, values_from = n_patients, values_fill = 0)
stable_pivot
# see which are the ID of these suspicious casese
sospetti_1_1 <- patients_stable %>%
  filter(stable_category == "1 -> 1", entry_research_group %in% c("EMCI", "MCI"))
sospetti_2_2 <- patients_stable %>%
  filter(stable_category == "2 -> 2", entry_research_group %in% c("CN", "SMC"))
sospetti_3_3 <- patients_stable %>%
  filter(stable_category == "3 -> 3", entry_research_group %in% c("CN", "EMCI", "MCI", "SMC", "LMCI"))
plot_sospetti <- function(data, sospetti, title) {
  data %>%
    filter(RID %in% sospetti$RID) %>%
    ggplot(aes(x = month, y = factor(RID), group = RID, color = factor(DIAGNOSIS))) +
    geom_line(alpha = 0.5) +
    geom_point(size = 1.5, alpha = 0.7) +
    facet_wrap(~entry_research_group, ncol = 1, scales = "free_y") + # separa per entry research group
    labs(
      title = title,
      x = "Month from baseline",
      y = "Patient ID",
      color = "Diagnosis"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
plot_1_1 <- plot_sospetti(data, sospetti_1_1, "Suspicious patients: 1 -> 1")
plot_2_2 <- plot_sospetti(data, sospetti_2_2, "Suspicious patients: 2 -> 2")
plot_3_3 <- plot_sospetti(data, sospetti_3_3, "Suspicious patients: 3 -> 3")
plot_1_1
plot_2_2
plot_3_3
#### Outlier analysis based on that ------------
order <- c("CN", "SMC", "EMCI", "LMCI", "MCI", "AD")
id <- sospetti_1_1$RID
id
sub <- data[data$RID %in% id, ]
sub$entry_research_group <- factor(sub$entry_research_group, levels = order)
sub=sub[order(sub$entry_research_group, sub$RID), ]
sub[,c("RID","entry","month","entry_research_group", "DIAGNOSIS")]
id <- sospetti_2_2$RID
id
sub <- data[data$RID %in% id, ]
sub$entry_research_group <- factor(sub$entry_research_group, levels = order)
sub=sub[order(sub$entry_research_group, sub$RID), ]
sub[,c("RID","entry","month","entry_research_group","DIAGNOSIS")]
id <- sospetti_3_3$RID
id
sub <- data[data$RID %in% id, ]
sub$entry_research_group <- factor(sub$entry_research_group, levels = order)
sub=sub[order(sub$entry_research_group, sub$RID), ]
sub[,c("RID","entry","nav","month","entry_research_group","DIAGNOSIS")]


# 0.5 Genotype distribution ------
temp <- data[, c("RID", "GENOTYPE")]
temp_unique <- unique(temp)
phase_count <- table(temp_unique$RID)
phase_count[phase_count > 1] #genotype non è variabile
#check that GENOTYPE doesn't change across patients: IT DOESN'T
categoria_per_paziente <- data %>%
  distinct(RID, entry_research_group, GENOTYPE)
table(categoria_per_paziente$entry_research_group, categoria_per_paziente$GENOTYPE)
ggplot(categoria_per_paziente, aes(x = GENOTYPE, fill = entry_research_group)) +
  geom_bar() +
  geom_text(
    aes(label = ..count..),
    stat = "count",
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  labs(
    title = "Distribution of genotypes by diagnostic group",
    x = "Genotype",
    y = "Number of patients",
    fill = "Diagnostic group"
  ) +
  theme_minimal(base_size = 14)
# /!\ should dig deeper to understand APOE <--> AD

# 0.6 MMSC distribution ------
groups <- c("CN","SMC","EMCI","LMCI","MCI","AD")  # ordina i gruppi come vuoi
for (g in groups) {
  data_subset <- data[data$entry_research_group == g, ]
  p <- ggplot(data_subset, aes(x = month, y = MMSCORE, group = RID, color = factor(RID))) +
    geom_line(alpha = 0.5) +
    geom_point(size = 1.5, alpha = 0.7) +
    labs(
      title = paste("MMSCORE trajectories -", g),
      x = "Month from baseline",
      y = "MMSE Score",
      color = "Patient ID"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")  

  print(p)
}
# /!\ might do also some more plots based on the diagnosis

# 0.7 PHASE distribution per patient and per research group ------
#overall
categoria_per_paziente <- data %>%
  distinct(RID, PHASE)
table(categoria_per_paziente$PHASE)
#per research group
temp <- data[, c("RID", "PHASE")]
temp_unique <- unique(temp)
phase_count <- table(temp_unique$RID)
phase_count[phase_count > 1] #phase è variabile
#check that PHASE doesn't change across patients: IT DOES
categoria_per_paziente <- data %>%
  group_by(RID, entry_research_group) %>%
  summarise(PHASE = max(PHASE), .groups = "drop") #take the last phase of the patient (can try also with min)
table(categoria_per_paziente$entry_research_group, categoria_per_paziente$PHASE)

# 0.8 Number of visits available distribution -------
library(dplyr)
library(ggplot2)

categorie_per_paziente <- data %>%
  group_by(RID) %>%
  summarise(n_categorie = n_distinct(entry), .groups = "drop")

ggplot(categorie_per_paziente, aes(x = factor(n_categorie))) +
  geom_bar(fill = "steelblue", color = "white") +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.3,
    size = 4
  ) +
  labs(
    title = "Distribution of Number of Different Entries Known per Patient",
    x = "Number of Different Entries",
    y = "Number of Patients"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(color = "black")
  )



# 0.9 Soumick Number of visits available distribution -------
pazienti_stabili <- data %>%
  group_by(RID) %>%
  filter(n_distinct(DIAGNOSIS) == 1) %>%
  ungroup()

visite_per_paziente <- pazienti_stabili %>%
  arrange(RID, month) %>%
  group_by(RID) %>%
  summarise(
    nav = n(),
    DIAGNOSIS = first(DIAGNOSIS),
    .groups = "drop"
  ) %>%
  mutate(
    DIAGNOSIS_label = case_when(
      DIAGNOSIS == 1 ~ "CN",
      DIAGNOSIS == 2 ~ "MCI",
      DIAGNOSIS == 3 ~ "AD"
    )
  )

totali <- visite_per_paziente %>%
  group_by(DIAGNOSIS_label) %>%
  summarise(n_totale = n(), .groups = "drop")

visite_per_paziente <- visite_per_paziente %>%
  left_join(totali, by = "DIAGNOSIS_label") %>%
  mutate(
    DIAGNOSIS_facet = paste0(DIAGNOSIS_label, " (N=", n_totale, ")")
  )

ggplot(visite_per_paziente, aes(x = factor(nav))) +
  geom_bar(fill = "steelblue") +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 4
  ) +
  facet_wrap(~ DIAGNOSIS_facet, ncol = 1) +
  labs(
    title = "Distribution of Number of Visits per Patient by Diagnosis (Stable Patients Only)",
    x = "Total number of visits",
    y = "Number of patients"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )



      


## follow-up length vs number of visits available ---------
followup_mnth <- data %>%
  group_by(RID) %>%
  summarise(
    n_visite = n_distinct(entry),
    durata_mesi = max(month),
    .groups = "drop"
  )
followup_mnth <- followup_mnth %>%
  left_join(data %>% distinct(RID, entry_research_group), by = "RID")
ggplot(followup_mnth, aes(x = durata_mesi, y = n_visite, color = entry_research_group)) +
  geom_count(alpha = 0.7) +
  labs(
    title = "Number of visits vs follow-up duration per patient",
    x = "Follow-up duration (months)",
    y = "Number of patients",
    color = "Diagnostic group",
    size = "Number of patients"
  ) +
  theme_minimal(base_size = 14)

## time lines per research group per patient ---------
data_grouped <- data %>%
  group_by(entry_research_group) %>%
  arrange(RID) %>%
  mutate(index = row_number()) %>%
  ungroup()
groups <- unique(data_grouped$entry_research_group)
for (g in groups) {
  data_subset <- data_grouped %>% filter(entry_research_group == g)
  
  p <- ggplot(data_subset, aes(x = month, y = factor(RID), group = RID, color = entry_research_group)) +
    geom_line(alpha = 0.5) +
    geom_point(size = 1, alpha = 0.7) +
    labs(
      title = paste("Timeline of visits -", g),
      x = "Month from baseline",
      y = "Patient ID",
      color = "Diagnostic group"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  
  print(p)
}

# 0.9 Age distribution -------
#overall
ggplot(data, aes(x = age)) +
  geom_histogram(
    binwidth = 5,
    fill = "steelblue",
    color = "white"
  ) +
  scale_x_continuous(
    breaks = seq(
      floor(min(data$age, na.rm = TRUE)),
      ceiling(max(data$age, na.rm = TRUE)),
      by = 5
    )
  ) +
  labs(
    title = "Age Distribution of All Patients",
    x = "Age (years)",
    y = "Number of patients"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 14),
    panel.grid.minor = element_blank()
  )


# per entry research group: it seems well stratified, no particular pattern
ggplot(data, aes(x = age, fill = entry_research_group)) +
  geom_histogram(binwidth = 5, color = "white", position="stack") +
  labs(
    title = "Age distribution by entry research group",
    x = "Age",
    y = "Number of patients",
    fill = "Entry research group"
  ) +
  theme_minimal(base_size = 14)
#time line per diagnosis
groups <- unique(data$entry_research_group)
for (g in groups) {
  data_subset <- data[data$entry_research_group == g, ]
  p <- ggplot(data_subset, aes(x = age, y = factor(RID), group = RID, color = factor(DIAGNOSIS))) +
    geom_line(alpha = 0.5) +
    geom_point(size = 1.5, alpha = 0.7) +
    labs(
      title = paste("Patient timelines by AGE -", g),
      x = "Age",
      y = "Patient ID",
      color = "Diagnosis"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
}

# Ideas that I want to try to implement 
#1.1 depth measures of each research group wrt CN and wrt AD. Compute risk = depth in cluster AD/depth in cluster CN
#1.2. box plots vs research group of the depth measures 
#1.3 scatterplot(mmse, risk)

























# Report PLOTS -------------
## Alluvional plot
library(dplyr)
library(ggplot2)
library(ggalluvial)

target_levels <- c("AD", "MCI", "CN")

plot_data <- data %>%
  group_by(RID) %>%
  arrange(age) %>% 
  summarise(
    Entry_cat = case_when(first(DIAGNOSIS) == 1 ~ "CN", 
                          first(DIAGNOSIS) == 2 ~ "MCI", 
                          first(DIAGNOSIS) == 4 ~ "AD"),
    Last_cat = case_when(last(DIAGNOSIS) == 1 ~ "CN", 
                         last(DIAGNOSIS) == 2 ~ "MCI", 
                         last(DIAGNOSIS) == 4 ~ "AD"),
    .groups = 'drop'
  ) %>%
  mutate(
    Entry_cat = factor(Entry_cat, levels = target_levels),
    Last_cat = factor(Last_cat, levels = target_levels)
  )

labels_entry <- plot_data %>%
  group_by(Entry_cat) %>%
  summarise(n_e = n()) %>%
  mutate(Entry_Label = paste0(Entry_cat, "\n(n = ", n_e, ")"))

labels_last <- plot_data %>%
  group_by(Last_cat) %>%
  summarise(n_l = n()) %>%
  mutate(Last_Label = paste0(Last_cat, "\n(n = ", n_l, ")"))

final_plot_data <- plot_data %>%
  left_join(labels_entry, by = "Entry_cat") %>%
  left_join(labels_last, by = "Last_cat")

final_plot_data$Entry_Label <- reorder(final_plot_data$Entry_Label, as.numeric(final_plot_data$Entry_cat))
final_plot_data$Last_Label <- reorder(final_plot_data$Last_Label, as.numeric(final_plot_data$Last_cat))

my_colors <- c("AD" = "#e31a1c", "CN" = "#33a02c", "MCI" = "#1f78b4")

ggplot(final_plot_data,
       aes(y = 1, axis1 = Entry_Label, axis2 = Last_Label)) +
  geom_alluvium(aes(fill = Entry_cat), width = 1/4, alpha = 0.45) +
  geom_stratum(width = 1/3, fill = "grey95", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5, fontface = "bold") +
  scale_fill_manual(values = my_colors) +
  scale_x_discrete(limits = c("Entry Visit", "Last follow-up"), expand = c(.1, .1)) +
  labs(title = "Diagnosis transitions from entry visit to last follow-up",
       y = "Number of Patients",
       fill = "Entry Diagnosis") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold")
  )
## Diagnostic path distr -------

data_clean <- data %>%
  filter(!is.na(DIAGNOSIS)) %>%
  mutate(
    DIAGNOSIS = case_when(
      DIAGNOSIS == 1 ~ "CN",
      DIAGNOSIS == 2 ~ "MCI",
      DIAGNOSIS == 4 ~ "AD",
      TRUE ~ as.character(DIAGNOSIS)
    )
  ) %>%
  arrange(RID, month)

stats_patients <- data_clean %>%
  group_by(RID) %>%
  summarise(
    path = paste(rle(as.character(DIAGNOSIS))$values, collapse = " -> "),
    n_states = length(rle(as.character(DIAGNOSIS))$values),
    .groups = "drop"
  ) %>%
  mutate(
    status_category = if_else(n_states == 1, "Stable", "Change")
  )

summary_table <- stats_patients %>%
  group_by(status_category, path) %>%
  summarise(patient_count = n(), .groups = "drop") %>%
  mutate(percentage = round((patient_count / sum(patient_count)) * 100, 2)) %>%
  arrange(desc(status_category), desc(patient_count))

total_patients <- nrow(stats_patients)
max_val <- max(summary_table$patient_count) * 1.3

cat("Total Number of Patients:", total_patients, "\n")
print(as.data.frame(summary_table))

ggplot(summary_table, aes(x = reorder(path, patient_count), y = patient_count, fill = status_category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(patient_count, " (", percentage, "%)")), 
            hjust = -0.1, 
            size = 3.5) +
  coord_flip(clip = "off") + 
  scale_y_continuous(limits = c(0, max_val)) +
  scale_fill_manual(values = c("Change" = "#E41A1C", "Stable" = "#377EB8")) +
  labs(
    title = "Diagnostic Path Distribution",
    subtitle = paste("Total N =", total_patients),
    x = "Diagnostic Path",
    y = "Number of Patients",
    fill = "Status"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 80, 10, 10)
  )




