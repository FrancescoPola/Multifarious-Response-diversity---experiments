rm(list = ls())
library(plotly)
library(tidyverse)
library(here)
library(patchwork)
library(DT)
library(mgcv)
library(gratia)
library(rlang)
library(vctrs)
library(scales)
library(broom)
library(reshape2)
library(ggtext)
library(ggsci)
library(ggpubr)
source.files <- list.files(here("r"), full.names = TRUE)
sapply(source.files, source, .GlobalEnv)

theme_set(theme_classic())

dd <- read.csv(here("data_microcosms/selection_mean.csv"))

# Define a function to map composition to species
map_composition_to_species <- function(composition) {
  species_mapping <- c("C" = "Colp", "L" = "Loxo", "D" = "Dexio",
                       "S" = "Spiro", "P" = "Para")
  species <- sapply(strsplit(composition, ""), function(x) {
    paste0(species_mapping[x], collapse = "_")
  })
  return(species)
}

# Apply the function to create the new column
dd <- dd %>%
  mutate(species = map_composition_to_species(composition))
# Remove "°C " from the temperature column
dd$temperature <- gsub(" °C ", "", dd$temperature)

# Remove " g/L" from the nutrients column
dd$nutrients <- gsub(" g/L", "", dd$nutrients)

df <- dd %>%
  slice(rep(seq_len(n()), each = 3)) %>%
  mutate(replicate = rep(1:3, times = nrow(dd))) %>%
  mutate(sample_ID = paste(composition, treatment, replicate, sep = "_"))

write.csv(df, here("data_microcosms/labels.csv"), quote = F, row.names = F)

# Group by nutrients and shuffle rows within each group
df_shuffled <- df %>%
  group_by(temperature) %>%
  slice(sample(n())) %>%
  ungroup()  # Remove grouping for the final dataframe
df_shuffled$temperature <- gsub("\\s*°C\\s*", "", df_shuffled$temperature)


write.csv(df_shuffled, here("data_microcosms/video.description.csv"), quote = F, row.names = F)
unique(df_shuffled$temperature)
