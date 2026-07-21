########################################################################################
########################################################################################
### Influenza Aerosol PFU/RNA Data 
### Joanna Pulit-Penaloza - A/Washington/2148/2025 (H5N5)
### Contextualizing H5N5 data with Nature Communications study
### 8 July 2026 - 20 July 2026
### Troy J. Kieran & Joanna Pulit-Penaloza
########################################################################################
########################################################################################
### Load packages
library(tidyverse)
library(tidylog)

########################################################################################
### Import Data

## Combined data from Nature Communications and Communications Biology JPP papers
data <- read.csv("combined_data_8July26.csv", header = TRUE, check.names = FALSE) 

########################################################################################

### Data Cleaning

## drop row number column
data <- data[, -1]

## remove the added NIOSH study viruses
## drop unneeded column
data <- data %>% 
  filter(sample_time != '2h_ cage') %>%
  select(-sample_time)

## rename viruses to add the /
data$Virus <- str_replace_all(data$Virus, "(\\D)(\\d)", "\\1/\\2")
data$Virus <- ifelse(data$Virus == "BC", "BC/PHL", data$Virus)

## set Sample & Transmission levels
data$Sample <- factor(data$Sample, levels = c("NW", "NIOSH", "SPOT"))
data$trans_cat <- factor(data$trans_cat, levels = c("non-transmissible", "transmissible"))
data$trans_3cat <- factor(data$trans_3cat, 
                          levels = c("Non-transmissible", "Low-transmissible", "High-transmissible"))

########################################################################################

### Calculate Values

data <- data %>%
  ## get unlogged values from log values
  mutate(PFU_unlog = if_else(is.na(PFU_unlog), 10^PFU, PFU_unlog),
         RNA_unlog = if_else(is.na(RNA_unlog), 10^RNA, RNA_unlog),
         day_num = case_when(day == "d1" ~ 1, day == "d2" ~ 2, day == "d3" ~ 3,
                             day == "d4" ~ 4, day == "d5" ~ 5)) %>% 
  ## get log values from unlogged values
  mutate(PFU = log10(PFU_unlog),
         RNA = log10(RNA_unlog)) %>%
  group_by(Virus, ferret, Sample) %>%
  ## calculate Area Under the Curve of unlogged values
  mutate(AUC_PFU = DescTools::AUC(day_num[day %in% c("d1", "d2", "d3")], 
                                  PFU_unlog[day %in% c("d1", "d2", "d3")], 
                                  method = "trapezoid", na.rm = TRUE),
         AUC_RNA = DescTools::AUC(day_num[day %in% c("d1", "d2", "d3")], 
                                  RNA_unlog[day %in% c("d1", "d2", "d3")], 
                                  method = "trapezoid", na.rm = TRUE),
         ## log the AUC values
         AUC_PFU = log10(AUC_PFU),
         AUC_RNA = log10(AUC_RNA)) %>%
  ungroup() %>%
  mutate(across(where(is.numeric), 
                ~ ifelse(is.nan(.) | is.na(.) | is.infinite(.), 0, .))) %>%
  filter(Sample != 'SPOT')

########################################################################################

### Transform to Long Format for Plotting

data_long <- data %>%
  group_by(Virus, ferret) %>%
  mutate(Sample2 = case_when(Sample == 'NW' ~ 'NW', 
                             Sample == 'NIOSH' ~ 'Air'),
         occurrence = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = Sample2, values_from = c(AUC_PFU, AUC_RNA),
              id_cols = c(Virus, Subtype, clade, genotype, day, ferret, trans_3cat)) %>%
  unnest(cols = everything()) %>% 
  distinct()

########################################################################################

### Make Ellipses Plots

## set group for ellipses
grp <- "trans_3cat" ## three (non/low/high) transmission categories

###

## consolidate non-2.3.4.4b genotypes to 'other'
data_long <- data_long %>%
  mutate(shape_grp = case_when(
    genotype %in% c("1", "H1", "H7", "H9") ~ "Other",
    TRUE ~ genotype))

###

### PFU plot

## calculate group mean for centroid points
cent_pfu <- data_long %>%
  group_by(.data[[grp]]) %>%
  summarise(mx = mean(AUC_PFU_NW), 
            my = mean(AUC_PFU_Air), .groups = "drop")

## drop perfectly collinear groups (i.e. H1)
data_ellip_pfu <- data_long %>%
  group_by(.data[[grp]]) %>%
  filter(abs(cor(AUC_PFU_NW, AUC_PFU_Air)) < 0.9999) %>%
  ungroup()

###

## plot group means with ellipses
ggplot() +
  ## individual points (optional)
  geom_point(data = data_long,
             aes(AUC_PFU_NW, AUC_PFU_Air, color = .data[[grp]], 
                 shape = shape_grp), alpha = 0.7, size = 3) +
  ## show spread of data with 95% CI
  stat_ellipse(data = data_ellip_pfu,
               aes(AUC_PFU_NW, AUC_PFU_Air, fill = .data[[grp]]),
               geom = "polygon", type = "norm", alpha = 0.15, level = 0.95) +
  ## mean centroid point - background highlight
  geom_point(data = cent_pfu,
             aes(mx, my), size = 4, color = 'red') +
  ## mean centroid point
  geom_point(data = cent_pfu,
             aes(mx, my, color = .data[[grp]]), size = 3) +
  scale_fill_viridis_d(end = 0.9) +
  scale_color_viridis_d(end = 0.9) +
  scale_shape_manual(values = c("D1.1" = 8, "B3.13" = 0, "A6" = 17)) +
  labs(x = "Nasal wash PFU AUC of days 1-3",
       y = "Airborne PFU AUC of days 1-3",
       color = NULL, fill = NULL, 
       shape = "2.3.4.4b Genotype") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text  = element_text(size = 14, face = "bold"),
        legend.position = "bottom",
        legend.position.inside = c(0.98, 0.001),
        legend.justification.inside = c(1, 0),
        legend.box = "vertical",
        legend.background = element_rect(fill = "white", color = NA)) +
  guides(fill = "none",
         color = guide_legend(
           title = "Transmissibility",
           position = "inside",
           ncol = 1,
           byrow = TRUE),
         shape = guide_legend(position = "bottom", nrow = 1)) -> plot_ellip_pfu

###

### RNA plot

## calculate group mean for centroid points
cent_rna <- data_long %>%
  group_by(.data[[grp]]) %>%
  summarise(mx = mean(AUC_RNA_NW), 
            my = mean(AUC_RNA_Air), .groups = "drop")

## drop perfectly collinear groups (i.e. H1)
data_ellip_rna <- data_long %>%
  group_by(.data[[grp]]) %>%
  filter(abs(cor(AUC_RNA_NW, AUC_RNA_Air)) < 0.9999) %>%
  ungroup()

###

ggplot() +
  ## individual points (optional)
  geom_point(data = data_long,
             aes(AUC_RNA_NW, AUC_RNA_Air, color = .data[[grp]], 
                 shape = shape_grp), alpha = 0.7, size = 3) +
  ## show spread of data with 95% CI
  stat_ellipse(data = data_ellip_rna,
               aes(AUC_RNA_NW, AUC_RNA_Air, fill = .data[[grp]]),
               geom = "polygon", type = "norm", alpha = 0.15, level = 0.95) +
  ## mean centroid point - background highlight
  geom_point(data = cent_rna,
             aes(mx, my), size = 4, color = 'red') +
  ## mean centroid point
  geom_point(data = cent_rna,
             aes(mx, my, color = .data[[grp]]), size = 3) +
  scale_fill_viridis_d(end = 0.9) +
  scale_color_viridis_d(end = 0.9) +
  scale_shape_manual(values = c("D1.1" = 8, "B3.13" = 0, "A6" = 17)) +
  labs(x = "Nasal wash RNA AUC of days 1-3",
       y = "Airborne RNA AUC of days 1-3",
       color = NULL, fill = NULL, 
       shape = "2.3.4.4b Genotype") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text  = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text  = element_text(size = 14, face = "bold"),
        legend.position = "bottom",
        legend.position.inside = c(0.98, 0.001),
        legend.justification.inside = c(1, 0),
        legend.box = "vertical",
        legend.background = element_rect(fill = "white", color = NA)) +
  guides(fill = "none",
         color = guide_legend(
           title = "Transmissibility",
           position = "inside",
           ncol = 1,
           byrow = TRUE),
         shape = guide_legend(position = "bottom", nrow = 1)) -> plot_ellip_rna


########################################################################################

### Save Plots to Image Files

## PFU
ragg::agg_tiff(filename = "plot_ellip_pfu_transmission_v4_JPP.tiff",
               width = 5.5, height = 5.5, units = 'in', res = 600,
               scaling = 1, compression = "lzw")
plot(plot_ellip_pfu)
invisible(dev.off())

## RNA
ragg::agg_tiff(filename = "plot_ellip_rna_transmission_v4_JPP.tiff",
               width = 5.5, height = 5.5, units = 'in', res = 600,
               scaling = 1, compression = "lzw")
plot(plot_ellip_rna)
invisible(dev.off())

##########################################################################################
##########################################################################################
### End of Code
##########################################################################################
##########################################################################################