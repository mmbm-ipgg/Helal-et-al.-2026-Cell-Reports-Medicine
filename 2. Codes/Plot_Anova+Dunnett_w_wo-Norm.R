#-----BC152 without normalization C0 and D0 + one-way ANOVA

# ---------- PRE-REQUISITE ---------------
#1 Chemin d'acces
path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))

#2 Creer une liste avec tous les fichiers .txt exportes par ImageJ
#file_list <- as.vector(t(read.table("FileList.txt", header=F, sep="\n")));
file_list = list.files(pattern = ".*.txt")

#3 Dataframe with all .txt files
dataset <- read.table(file_list, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
dataset <- as.data.frame(dataset)

# ---------- PACKAGES ---------------
#Packages
library(ggplot2)
library(scales)
library(cowplot)
library(MASS)
library(plyr)
library(dplyr)
library(sm)
library(beeswarm)
library(ggsignif)
library(stringr)
library(ggpubr)
library(rstatix)
library(ggpattern)

# ---------- DATA ORGANIZATION ---------------
#1 Mean data for 9 values of live/dead
dataset <- dataset %>%
  group_by(Drug, Tumor, Replica, Day, Concentration) %>% 
  summarise("mean_dead" = mean(Number_dead), "mean_live" = mean(Number_live))

#2 Cell viability (live/total)
dataset = cbind(dataset, "viability" = dataset$mean_live/(dataset$mean_dead+dataset$mean_live))

#4 Mean of the technical replicates
dataset_mean_bio <- dataset %>%
  group_by(Drug, Tumor, Day, Concentration) %>% 
  summarise("mean_viability" = mean(viability), "sd_viability" = sd(viability))

#5 Mean of the biological replicates (for y position of mean)
dataset_mean <- dataset_mean_bio %>%
  group_by(Drug, Day, Concentration) %>% 
  summarise("mean_viability_bio" = mean(mean_viability), "sd_viability" = sd(mean_viability))


#--NORMALIZED BY C0--
#3bis Relative cell viability normalized by the mean of C0 for each day from 0 to day 4 (Which gives that Control at D0 = 1)
#For loop to normalized all data at different days over the concentration = 0
results_list <- list()
for (d in 0:4) {
  # Filtrer les données pour le jour actuel
  day <- dataset[dataset$Day == as.character(d), ]
  # Calculer les valeurs normalisées
  day_C0 <- cbind(
    day,
    "viability_C0" = day$viability / mean(day$viability[day$Drug == "Control"])
  )
  # Sauvegarder dans la liste
  results_list[[paste0("Day", d)]] <- day_C0
}

#Combiner tous les jours dans un seul data.frame
dataset_normalizedC0 <- do.call(rbind, results_list)

#4bis Mean of the technical replicates (donc tous les points sur le grpahes seront des replicas biologiques)
dataset_mean_bio_normalizedC0 <- dataset_normalizedC0 %>%
  group_by(Drug, Tumor, Day, Concentration) %>% 
  summarise("mean_viability_normalizedC0" = mean(viability_C0), "sd_viability_normalizedC0" = sd(viability_C0))

#5 Mean of the biological replicates (for y position of mean)
dataset_mean_normalizedC0 <- dataset_mean_bio_normalizedC0 %>%
  group_by(Drug, Day, Concentration) %>% 
  summarise("mean_viability_bio_normalizedC0" = mean(mean_viability_normalizedC0), "sd_viability_normalizedC0" = sd(mean_viability_normalizedC0))

# ---------- STATISTIQUES WILCOXON ---------------
#1 Wilcoxon on all points
stat.test_wilcoxon <- compare_means(
  formula = mean_viability ~ Concentration,
  data = dataset_mean_bio,
  group.by = "Day",
  method = "wilcox.test",
  paired=F
)

#2 Keep only the comparison with DAY 0
stat.test_wilcoxon = stat.test_wilcoxon[stat.test_wilcoxon$group2=='0', ]

#3 Find x and y position for further plotting stats in front of the correct dataset
stat.test_wilcoxon <- stat.test_wilcoxon %>%
  left_join(
    dataset_mean %>%
      mutate(Concentration = as.character(Concentration)) %>%
      select(Day, Concentration, mean_viability_bio),
    by = c("Day", "group1" = "Concentration")
  ) %>%
  mutate(
    x_position = Day + 0.5,
    y_position = mean_viability_bio
  )

#4 Change the name of group1 in "concentration"
stat.test_wilcoxon <- stat.test_wilcoxon %>%
  rename(Concentration = group1)

#--NORMALIZED BY C0--
#1bis Wilcoxon on all datasets
stat.test_wilcoxon_normalizedC0 <- compare_means(
  formula = mean_viability_normalizedC0 ~ Concentration,
  data = dataset_mean_bio_normalizedC0,
  group.by = "Day",
  method = "wilcox.test",
  paired=F
)

#2bis Keep only the comparison with DAY 0
stat.test_wilcoxon_normalizedC0 = stat.test_wilcoxon_normalizedC0[stat.test_wilcoxon_normalizedC0$group2=='0', ]

#3bis Find x and y position for further plotting stats in front of the correct dataset
stat.test_wilcoxon_normalizedC0 <- stat.test_wilcoxon_normalizedC0 %>%
  left_join(
    dataset_mean_normalizedC0 %>%
      mutate(Concentration = as.character(Concentration)) %>%
      select(Day, Concentration, mean_viability_bio_normalizedC0),
    by = c("Day", "group1" = "Concentration")
  ) %>%
  mutate(
    x_position = Day + 0.5,
    y_position = mean_viability_bio_normalizedC0
  )

#4bis Change the name of group1 in "concentration"
stat.test_wilcoxon_normalizedC0 <- stat.test_wilcoxon_normalizedC0 %>%
  rename(Concentration = group1)

# ---------- STATISTICS ONE-WAY ANOVA ---------------
#1 Remove the Day0 as it contains only one condition (C0) and anova compares several conditions for each day 
dataset_woD0 <- dataset_mean_bio %>%
  group_by(Day) %>%
  filter(n_distinct(Concentration) >= 2) %>%
  ungroup()

#2 Check of SD are different
library(car)
Levene_test <- dataset_woD0 %>%
  group_by(Day) %>%
  group_modify(~ {
    # Ensure Concentration is treated as a factor
    .x$Concentration <- as.factor(.x$Concentration)
    
    # Run Levene's test
    test <- LeveneTest(mean_viability ~ Concentration, data = .x)
    
    # Convert to tibble to return
    as_tibble(test)
  }) %>%
  ungroup()
Levene_test

#2 Anova on all the datapoints
stat.test_anova <- compare_means(
  formula = mean_viability ~ Concentration,
  data = dataset_woD0,
  group.by = "Day",
  method = "anova"
)

#3 IF anova is significant, run the Post-hoc test (ANOVA tells you that there is a difference, not where.)
library(dplyr)
library(DescTools)

#add the significant legend
signif_codes <- function(p) {
  cut(
    p,
    breaks = c(-Inf, 0, 0.001, 0.01, 0.05, 0.1, 1),
    labels = c("****", "***", "**", "*", "ns", "ns"),
    right = TRUE
  )
}

posthoc_Dunnett <- dataset_woD0 %>%
  group_by(Day) %>%
  group_modify(~ {
    dt <- DunnettTest(
      x = .x$mean_viability,
      g = .x$Concentration,
      control = "0"
    )
    
    out <- as.data.frame(dt[[1]])
    out$Comparison <- rownames(out)
    
    # Add Significance codes
    out$Signif <- signif_codes(out$pval)
    
    out
  }) %>%
  ungroup()
posthoc_Dunnett

#4 Trouver les x et y positions pour annoter les p-values à coté des points de D1-D4
posthoc_Dunnett <- posthoc_Dunnett %>%
  # Extract treatment from Comparison
  mutate(Treatment = sub("-0$", "", Comparison)) %>%
  # Join the mean viability
  left_join(
    dataset_mean %>%
      mutate(Concentration = as.character(Concentration)),
    by = c("Day", "Treatment" = "Concentration")
  ) %>%
  # Add positions for plotting
  mutate(
    x_position = Day + 0.5,
    y_position = mean_viability_bio
  )

#5 Changer le nom de de la colonne Treatment en Concentration pour avoir les même couleurs des stat et des points dans ggplot
posthoc_Dunnett <- posthoc_Dunnett %>%
  rename(Concentration = Treatment, p.signif = Signif)


#--NORMALIZED BY C0--
#1bis Remove the Day0 as it contains only one condition (C0) and anova compares several conditions for each day 
dataset_woD0_normalizedC0 <- dataset_mean_bio_normalizedC0 %>%
  group_by(Day) %>%
  filter(n_distinct(Concentration) >= 2) %>%
  ungroup()

#2 Anova sur tous les time points
stat.test_anova_normalizedC0 <- compare_means(
  formula = mean_viability_normalizedC0 ~ Concentration,
  data = dataset_woD0_normalizedC0,
  group.by = "Day",
  method = "anova"
)

#3 IF anova is significant, run the Post-hoc test (ANOVA tells you that there is a difference, not where.)
library(dplyr)
library(DescTools)
#add the significant legend
signif_codes <- function(p) {
  cut(
    p,
    breaks = c(-Inf, 0, 0.001, 0.01, 0.05, 0.1, 1),
    labels = c("****", "***", "**", "*", "ns", "ns"),
    right = TRUE
  )
}

posthoc_Dunnett_normalizedC0 <- dataset_woD0_normalizedC0 %>%
  group_by(Day) %>%
  group_modify(~ {
    dt <- DunnettTest(
      x = .x$mean_viability_normalizedC0,
      g = .x$Concentration,
      control = "0"
    )
    
    out <- as.data.frame(dt[[1]])
    out$Comparison <- rownames(out)
    
    # Add Significance codes
    out$Signif <- signif_codes(out$pval)
    
    out
  }) %>%
  ungroup()
posthoc_Dunnett_normalizedC0

#4 Find x and y position for further plotting stats in front of the correct datapoints
posthoc_Dunnett_normalizedC0 <- posthoc_Dunnett_normalizedC0 %>%
  # Extract treatment from Comparison
  mutate(Treatment = sub("-0$", "", Comparison)) %>%
  # Join the mean viability
  left_join(
    dataset_mean_normalizedC0 %>%
      mutate(Concentration = as.character(Concentration)),
    by = c("Day", "Treatment" = "Concentration")
  ) %>%
  # Add positions for plotting
  mutate(
    x_position = Day + 0.5,
    y_position = mean_viability_bio_normalizedC0
  )

#5 CHange the column name "Treatment" in "Concentration"
posthoc_Dunnett_normalizedC0 <- posthoc_Dunnett_normalizedC0 %>%
  rename(Concentration = Treatment, p.signif = Signif)

# ---------- ESTHETICS ---------------
#1. Theme
theme_Lea <- theme(
  plot.title = element_text(size = 18, face = "bold", hjust = 0.5, lineheight = 1.2),   # Title font size, bold, centered, line spacing
  panel.grid.major = element_blank(),                                                   # Remove major gridlines
  panel.grid.minor = element_blank(),                                                   # Remove minor gridlines
  panel.border = element_blank(),                                                       # Remove panel border
  panel.background = element_rect(fill = 'transparent'),                                 # Transparent panel background
  axis.text.x = element_text(size = 15, color = 'black'),                                # X-axis tick labels
  axis.text.y = element_text(size = 15, color = 'black'),                                # Y-axis tick labels
  axis.title.x = element_text(size = 15),                                               # X-axis title
  axis.title.y = element_text(size = 15),                                               # Y-axis title
  axis.line = element_line(color = "black", linewidth = 0.6),                            # Axis lines: black, thickness 0.6
  axis.ticks = element_line(color = "black", linewidth = 0.6),                           # Tick marks: black, same thickness as axis
  legend.position = "right",                                                            # Legend position on the right
  legend.title = element_blank(),                                                       # Remove legend title
  legend.key = element_blank(),                                                         # Remove legend key background
  legend.text = element_text(size = 12),                                               # Legend text size
  strip.background = element_rect(fill = "white", color = "black", linewidth = 1),       # Facet strip: white background, black border
  strip.text = element_text(size = 16, face = "bold"),                                   # Facet label text: bold, size 16
  strip.placement = "outside",                                                          # Place facet strips outside plot
  strip.text.x = element_text(margin = margin(t = 3, b = 3))                              # Facet strip text margin top/bottom
)

#2 Colors
library(RColorBrewer)
base_color <- "#08306B"
gradient_colors <- colorRampPalette(c("white", base_color))(5)
my_colors = tail(gradient_colors, 5)

# ---------- PLOTS ---------------
#1 Cell viability without any normalization
ggplot(data = dataset_mean_bio, aes(y = mean_viability*100, x = Day,
                                    color = ifelse(Drug == "Control", "Control", as.character(Concentration)))) +
  geom_point(size = 2, alpha = 0.3) +                                                                                # raw points
  stat_summary(fun.data = mean_se, geom = "linerange", size = 0.5) +                                # SE bars
  stat_summary(fun = "mean", geom = "point", size = 4) +                                          # colored mean
  geom_text(data = posthoc_Dunnett,
            aes(x = x_position, y = y_position*100, label = p.signif, color = as.factor(Concentration)),
            size=5,
            show.legend = F) +
  scale_color_manual(
    values = c("Control" = "grey70",
               setNames(my_colors, sort(unique(dataset_mean_bio$Concentration))))) +
  labs(x = "Time (days)", y = "Cell viability (%)", color = "[Drug] (µM)") +
  ggtitle("5FU") +
  theme_Lea +
  theme(legend.title=element_text(size=12, color='black')) +
  ylim(c(0, 150))
ggsave("Viability_5FU_C0_D0_PHDunnett.pdf", width = 14, height = 10, units = "cm")
ggsave("Viability_5FU_C0_D0_PHDunnett.png", width = 14, height = 10, units = "cm")

#2 Cell viability without any normalization for the figure S2
dataset_mean_bio = dataset_mean_bio[!dataset_mean_bio$Concentration==1 & !dataset_mean_bio$Concentration==20, ]
ggplot(data = dataset_mean_bio, aes(y = mean_viability*100, x = Day,
                                    color = ifelse(Drug == "Control", "Control", as.character(Concentration)))) +
  geom_point(size = 2, alpha = 0.3) +                                                                                # raw points
  stat_summary(fun.data = mean_se, geom = "linerange", size = 0.5) +                                # SE bars
  stat_summary(fun = "mean", geom = "point", size = 4) +                                          # colored mean
  scale_color_manual(
    values = c("Control" = "grey70", "10" = "#6baed6cc", "50" = "#08306bff")) +
  labs(x = "Time (days)", y = "Cell viability (%)", color = "[Drug] (µM)") +
  ggtitle("5FU") +
  theme_Lea +
  scale_y_continuous(
    limits = c(0, 125),
    breaks = seq(0, 125, by = 25)
  ) +
  theme(legend.title=element_text(size=12, color='black'))
ggsave("Viability_5FU_C0_D0_10-50.pdf", width = 14, height = 10, units = "cm")
ggsave("Viability_5FU_C0_D0_10-50.png", width = 14, height = 10, units = "cm")

#3 Relative cell viability normalized by the mean C0 of each day (Which gives that Control at D0 = 1)
ggplot(data = dataset_mean_bio_normalizedC0, aes(y = mean_viability_normalizedC0, x = Day,
                                    color = ifelse(Drug == "Control", "Control", as.character(Concentration)))) +
  geom_point(size = 2, alpha = 0.3) +                                                              # raw points
  stat_summary(fun.data = mean_se, geom = "linerange", size = 0.5) +                                # SE bars
  stat_summary(fun = "mean", geom = "point", size = 4) +                                          # colored mean
  geom_text(data = posthoc_Dunnett_normalizedC0,
            aes(x = x_position, y = y_position, label = p.signif, color = as.factor(Concentration)),
            size=5,
            show.legend = F) +
  scale_color_manual(
    values = c("Control" = "grey70",
               setNames(my_colors, sort(unique(dataset_mean_bio_normalizedC0$Concentration))))) +
  labs(x = "Time (days)", y = "Relative cell viability\n(Normalized by control)", color = "[Drug] (µM)") +
  ggtitle("5FU") +
  theme_Lea +
  theme(legend.title=element_text(size=12, color='black')) +
  ylim(c(0, 1.5))
ggsave("Viability_5FU_D0_normbyC0_PHDunnett.pdf", width = 14, height = 10, units = "cm")
ggsave("Viability_5FU_D0_normbyC0_PHDunnett.png", width = 14, height = 10, units = "cm")

# ---------- SUMMARY ---------------
#1 Mean of the data of dataset_mean_bio_normalizedC0 per day and per concentration
summary_df <- dataset_mean_bio %>%
  group_by(Day, Concentration) %>%
  summarise(
    mean_mean_viability = mean(mean_viability, na.rm = TRUE),
    sd_mean_viability   = sd(mean_viability, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
summary_df

#2 Counts of all the tumor per model
tumor_per_model <- dataset_mean_bio %>%
  group_by(Drug) %>%
  summarise(
    n_tumors = n_distinct(Tumor),
    .groups = "drop"
  )
N=sum(tumor_per_model$n_tumors[tumor_per_model$Drug=="5FU"])

#3 Counts of all the gels (technical replicates) per model
number_gels <- dataset %>%
  group_by(Drug) %>%
  summarise(
    n_gels = n_distinct(Replica, Tumor),
    .groups = "drop"
  )
n=sum(number_gels$n_gels)