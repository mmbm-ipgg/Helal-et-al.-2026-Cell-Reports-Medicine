#-----Biopsy human - Control and drug

# ---------- PRE-REQUISITE ---------------
#Chemin d'acces
path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))

#Create a list with all datasets
file_list = list.files(pattern = ".*.txt")

dataset <- read.table(file_list, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
dataset <- as.data.frame(dataset)

# ---------- PACKAGES ---------------
library(ggplot2)
library(scales)
library(cowplot)
library(MASS)
library(plyr)
library(sm)
library(beeswarm)
library(ggsignif)
library(dplyr)
library(stringr)
library(ggpubr)
library(rstatix)
library(ggpattern)
library(grDevices)
library(RColorBrewer)

# ---------- DATA ORGANIZATION ---------------
#1 Mean data for 9 values of live/dead
dataset <- dataset %>%
  group_by(Drug, Tumor, Replica, Day, Concentration) %>% 
  summarise("mean_dead" = mean(N_Dead), "mean_live" = mean(N_Live))

#2 Cell viability (live/total)
dataset = cbind(dataset, "viability" = dataset$mean_live/(dataset$mean_dead+dataset$mean_live))

#3 Remove the tumor that are <50% at D0
dataset <- dataset %>%
  group_by(Tumor) %>%
  filter(
    !any(Day == 0 & viability < 0.5, na.rm = TRUE)
  ) %>%
  ungroup()

#4 Reorder the number (name) of the tumors
dataset <- dataset %>%
mutate(
  Tumor_reorder = factor(
    dplyr::recode(
      as.character(Tumor),
      "2"  = "1",
      "3"  = "7 (Biopsy)",
      "4"  = "2",
      "6"  = "3",
      "7"  = "4",
      "8"  = "5",
      "13" = "6"
    ),
    levels = c("1", "2", "3", "4", "5", "6", "7 (Biopsy)")
  )
)

#5 Only D4, or worst case D3
dataset_D4 <- dataset[!dataset$Day == "0", ]
dataset_D4 <- dataset_D4[!dataset_D4$Day == "1", ]

#6 Only D0
dataset_D0 <- dataset[dataset$Day == "0", ]

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

#2. colors
base_color <- "#762A83"
gradient_colors <- colorRampPalette(c("#FD8D3C", base_color))(7)
my_colors = tail(gradient_colors, 7)
display.brewer.all()
brewer.pal(n = 9, name = "Oranges")


# ---------- PLOTS ---------------
#1 Cell viability drugs vs. control at D4
ggplot(data = dataset_D4, aes(y = viability*100, x = Drug,
                                     color = ifelse(Drug == "Control", "Control", as.character(Tumor_reorder)))) +
  geom_point(size=4.5, alpha=0.4) +
  stat_summary(fun.data = mean_se, geom = "linerange", size = 0.5, color = "black") +                                # SE bars
  stat_summary(fun = "mean", geom = "point", size = 4, color = "black") +                                            # big black mean
  stat_summary(fun = "mean", geom = "point", size = 3) +                                                             # colored mean
  labs(x = "", y = "Cell viability (Day 4) (%)", color = "# Tumor") +
  theme_Lea +                                                                                                       = 0.5) +
  theme(legend.title=element_text(size=12, color='black')) +
  scale_x_discrete(
    limits = function(x) c("Control", setdiff(x, "Control"))
  ) +
  scale_color_manual(
    values = c(
      "Control" = "grey70",
      setNames(my_colors, sort(unique(dataset_D4$Tumor_reorder))))) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1),   # rotate labels
  ) +
  facet_grid(~ Tumor_reorder, scales = "free_x", space = "free_x") +
  ylim(0, 100)
ggsave("Human_biopsy_C0.pdf", width = 28, height = 10, units = "cm")
ggsave("Human_biopsy_C0.png", width = 28, height = 10, units = "cm")


#2 Cell viability of the controls at D0
ggplot(data = dataset_D0, aes(y = viability*100, x = Drug, color=as.factor(Tumor_reorder), group=interaction(Tumor_reorder, Drug))) +
  geom_point(size=4.5, alpha=0.4) +
  stat_summary(fun.data = mean_se, geom = "linerange", size = 0.5, color = "black") +                                # SE bars
  stat_summary(fun = "mean", geom = "point", size = 4, color = "black") +                                            # big black mean
  stat_summary(fun = "mean", geom = "point", size = 3) +                                                             # colored mean
  labs(x = "", y = "Cell viability (Day 0) (%)", color = "# Tumor") +
  theme_Lea +
  geom_hline(yintercept = 50, linetype = "dashed", color = "black", size = 0.5) +
  theme(legend.title=element_text(size=12, color='black')) +
  scale_color_manual(values = c(my_colors)) +
  ylim(0, 100)
ggsave("Human_biopsy_Control_D0.pdf", width = 8, height = 10, units = "cm")
ggsave("Human_biopsy_Control_D0.png", width = 8, height = 10, units = "cm")