########################################################################################
########################################################################################
## Troy J. Kieran
## 10 March 2026 - 1 June 2026
## Joanna Pulit-Penaloza - A/Washington/2148/2025 (H5N5)
## Ferret Symptom Prevalence Figure
########################################################################################
########################################################################################

## Load Libraries
library(tidyverse)
library(tidylog)
#library(extrafont) ## not needed as Arial seems like default

########################################################################################

## Import and load fonts
#font_import()
#loadfonts(device = "win")

########################################################################################

## Import data
data <- read.csv("H5N5_symptom_prevalence_10Mar2026.csv")

########################################################################################

## dataframe for lethal annotation bar & labels
deaths <- data.frame(
  Day = c(3, 7),   ## days when individuals died
  count = c(1, 2)) ## number that died

## heatmap of symptom prevalence
data %>%
  group_by(Virus, Day, symptom) %>%
  reframe(sum = sum(symptom_present, na.rm = TRUE)) %>%
  mutate(sum = factor(sum, levels = c(0, 1, 2, 3))) %>%
  ggplot(aes(x = factor(Day), y = symptom, fill = sum)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_d(option = "viridis", end = 0.9, 
                       direction = -1, drop = TRUE) +
  labs(x = "Day Post-inoculation", y = NULL, 
       fill = "Number of Ferrets") +
  scale_y_discrete(labels = c(
    "lethargy" = "Lethargy",
    "lack of grooming" = "Lack of Grooming",
    "grimace" = "Grimace",
    "Nas Dis" = "Nasal Discharge",
    "snz" = "Sneeze",
    "hunched" = "Hunched",
    "diah" = "Diarrhea")) +
  theme_void() +
  theme(strip.background = element_rect(fill = "white"),
        text = element_text(family = "Arial"),
        strip.text = element_text(face = "bold"),
        legend.position = 'bottom',
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold")) +
  ## when legend was contiuous color, non-factor sum
  # guides(fill = guide_colorbar(
  #   barwidth  = unit(5, "cm"),  ## make longer, when at bottom
  #   barheight = unit(0.2, "cm"), ## make taller, when at bottom
  #   ticks = TRUE)) +
  guides(fill = guide_legend(nrow = 1)) +
  ## add the annotation bar for lethal
  geom_tile(data = deaths, 
            aes(x = factor(Day), y = 8, height = 1, width = 1), 
            fill = "red3", alpha = 0.7) +
  ## add numbers for lethal count to annotation bar
  geom_text(data = deaths, 
            aes(x = factor(Day), y = 8, label = count),
            color = "white", size = 5, fontface = "bold", 
            vjust = 0.5, inherit.aes = FALSE) #-> plot_symptom_prevalence

# ragg::agg_tiff(filename = "plot_symptom_prevalence.tiff",
#                width = 6, height = 3, units = 'in', res = 600,
#                scaling = 0.9, compression = "lzw")
# plot(plot_symptom_prevalence)
# invisible(dev.off())

########################################################################################
########################################################################################
### End of Code
########################################################################################
########################################################################################
