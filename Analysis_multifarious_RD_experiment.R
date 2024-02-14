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
library(mgcv)
source.files <- list.files(here("r"), full.names = TRUE)
sapply(source.files, source, .GlobalEnv)



# In this code, we are going to fit linear regressions to calculate 
# carry capacity (K) and growth rate (r) for each species and treatment combination


dd <- read.csv(here("data/dd_reaction_norms.csv"))


# Convert Date to Date format (if not already)
dd$date <- as.Date(dd$date)

# Arrange the data frame by Species and Date
dd <- dd %>% arrange(species, date)


# Create a reference data frame with all combinations of Species and Date
reference_data <- expand.grid(species = unique(dd$species), Date = unique(dd$date))

# Merge the reference data with your actual data
dd <- merge(reference_data, dd, all.x = TRUE)

# Add a new column 'day' indicating the number of the sampling for each species
dd <- dd %>%
  group_by(species) %>%
  mutate(day = match(date, unique(date))) %>% 
  select(-Date)

write.csv(dd, file = here("data/data_reaction_norms.csv"), quote = FALSE, row.names = FALSE)


dd <- read_csv(here("data/data_reaction_norms.csv"))


# convert nutrients to factor
dd$nutrients <- as.factor(dd$nutrients)
# plot sp densities over time

ggplot(dd, aes(x = day, y = log(mean.dens.ml), color = nutrients)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +  # You can adjust the smoothing method and add confidence intervals if needed
  scale_color_viridis_d(option = "vulcano") +  # Using viridis vulcano color palette
  facet_grid(temperature ~ species) +
  theme_bw() +
  scale_x_continuous(breaks = seq(min(dd$day), max(dd$day), by = 1)) 


## try with one species
# look at exponential growth phase for one temperature
p_dexio_28 <- dd %>% 
  filter(species == "Dexiostoma", day <= 5, temperature == 28) %>% 
  ggplot(aes(x = day, y = log(mean.dens.ml + 0.001), color = nutrients)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +  # You can adjust the smoothing method and add confidence intervals if needed
  scale_color_viridis_d(option = "magma") +
  theme_bw() +
  scale_x_continuous(breaks = seq(min(dd$day), max(dd$day), by = 1))
  

dexio_28 <- dd  %>% 
  filter(species == "Dexiostoma", day <= 5, temperature == 28) %>% 
  rename(density = "mean.dens.ml")


dexio_28 <- dexio_28 %>% group_by(nutrients, sample_ID)

### Calculate carrying capacity (K) as the highest population density for each species and treatment 

carrying_capacity_dexio_28<- dexio_28 %>%
  summarise(
    CarryingCapacity = max(density),
    .groups = 'drop'
  )

### Calculate growth rate as the slope of the linear regression of density ~ time for each species and treatment combination

# Add a small constant to avoid issues with log(0)
dexio_28$density <- dexio_28$density + 0.001

# Calculate growth rates for each nutrient level
growth_rates_dexio_28 <- dexio_28 %>%
  summarise(
    Intercept = coef(lm(log(density) ~ day))[1],
    GrowthRate = coef(lm(log(density) ~ day))[2],
    RSquared = summary(lm(log(density) ~ day))$r.squared
  ) %>%
  ungroup()  # Drop the grouping information

dexio_28_r_K <- merge(growth_rates_dexio_28, carrying_capacity_dexio_28, by = c( "nutrients", "sample_ID"))

head(dexio_28_r_K)
dexio_28_r_K$nutrients <- as.integer(dexio_28_r_K$nutrients)
# Fit a GAM for growth rate
gam_model <- gam(GrowthRate ~ s(nutrients, k = 5),
                 data = dexio_28_r_K,
                 method = "REML")


# Display the summary of the GAM model
summary(gam_model)



new_data <- expand_grid(nutrients = seq(1, 5, by= 0.02))

predicted <- predict(gam_model, newdata = new_data)

rates <- cbind(new_data, predicted)


p1 <- draw(gam_model, rug = FALSE, dist = 0.5)
p_dexio_28 + p1





# Dexio all Ts 

dexio_p <- dd %>% 
  filter(species == "Dexiostoma", day <= 5,) %>% 
  ggplot(aes(x = day, y = log(mean.dens.ml + 1), color = as.factor(nutrients))) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +  # You can adjust the smoothing method and add confidence intervals if needed
  scale_color_viridis_d(option = "magma") +
  facet_grid( ~ temperature ) +
  theme_bw() +
  scale_x_continuous(breaks = seq(min(dd$day), max(dd$day), by = 1))

## for temperature <= 24 not linear relationship...
dexio <- dd %>% filter(species == "Dexiostoma")
# Group the data by species, nutrient level, temperature, and sample ID
grouped_dexio <- dexio %>%
  group_by(nutrients, temperature, sample_ID) %>% 
  rename(density = mean.dens.ml)

### Calculate carrying capacity (K) as the highest population density for each species and treatment 

carrying_capacity_dexio <- grouped_dexio %>%
  summarise(
    CarryingCapacity = max(density),
    .groups = 'drop'
  )

### Calculate growth rate as the slope of the linear regression of density ~ time for each species and treatment combination

grouped_dexio_d6 <- dexio %>% filter(day <= 6) %>% 
  group_by(nutrients, temperature, sample_ID) %>% 
  rename(density = mean.dens.ml)

growth_rates_grouped_dexio <- grouped_dexio_d6 %>%
  summarise(
    Intercept = coef(lm((density + 0.001) ~ day))[1],
    GrowthRate = coef(lm((density+ 0.001 ) ~ day))[2],
    RSquared = summary(lm((density + 0.001) ~ day))$r.squared,
    .groups = 'drop'
  )

# Display the resulting carrying capacity and growth rate data frames
grouped_dexio_r_K <- merge(growth_rates_grouped_dexio, carrying_capacity_dexio, by = c("nutrients", "temperature", "sample_ID"))
head(grouped_dexio_r_K)

grouped_dexio_r_K$nutrients <- as.integer(grouped_dexio_r_K$nutrients)
# Fit a GAM for growth rate
gam_model_dexio_grouped <- gam(GrowthRate ~ti(temperature) + ti(nutrients) + te(temperature, nutrients, k = c(5, 5)),
                 data = grouped_dexio_r_K, method = "REML")


# Display the summary of the GAM model
summary(gam_model_dexio_grouped)



new_data <- expand_grid(temperature = seq(18, 28, by= 0.02),
                        nutrients = seq(1, 5, by= 0.01))

predicted <- predict(gam_model_dexio_grouped, newdata = new_data)

rates <- cbind(new_data, predicted)
rates <- rates %>%
  relocate(temperature, nutrients, predicted)
head(rates)

# plot surface
p2 <- draw(gam_model_dexio_grouped, rug = FALSE, dist = 0.5)


dexio_p + p2





# Now we apply the code same code to all species
dd <- read.csv(here("Data/data_reaction_norms.csv"))
head(dd)
unique(dd$day)

dd <- dd %>% 
  rename(density = "mean.dens.ml")

# create a data frame with only the days of exponential growth (first 6 days)
dd_d6 <- dd %>% filter(day <= 6) 


# Group the data by species, nutrient level, temperature, and sample ID
grouped_data <- dd %>%
  group_by(species, nutrients, temperature, sample_ID)

### Calculate carrying capacity (K) as the highest population density for each species and treatment 

carrying_capacity_df <- grouped_data %>%
  summarise(
    CarryingCapacity = max(density),
    .groups = 'drop'
  )

### Calculate growth rate as the slope of the linear regression of density ~ time for each species and treatment combination

grouped_data_d6 <- dd_d6 %>%
  group_by(species, nutrients, temperature, sample_ID)

growth_rates_df <- grouped_data_d6 %>%
  summarise(
    Intercept = coef(lm(log(density + 0.001) ~ day))[1],
    GrowthRate = coef(lm(log(density + 0.001) ~ day))[2],
    RSquared = summary(lm(log(density + 0.001) ~ day))$r.squared,
    .groups = 'drop'
  )

# Display the resulting carrying capacity and growth rate data frames
dd_r_K <- merge(growth_rates_df, carrying_capacity_df, by = c("species", "nutrients", "temperature", "sample_ID"))
head(growth_rates_df)




# List to store ggplot objects
ggplot_objects <- list()

# Loop through each species
for (species_name in unique(dd$species)) {
  
  # Filter data for the current species
  subset_data <- dd %>%
    filter(species == species_name)
  
  # Create ggplot object
  ggplot_object <- ggplot(subset_data, aes(x = day, y = log(density + 1), color = as.factor(nutrients))) +
    geom_point(data = subset_data, aes(x = day, y = log(density + 1), color = as.factor(nutrients))) +
    geom_smooth(data = subset_data, method = "loess", se = FALSE) +
    geom_smooth(data = subset_data %>% filter(day <= 5), method = "lm", se = FALSE, linetype = "dashed", color = "black", aes(group = nutrients, color = nutrients)) +  # Add LM line
    scale_color_viridis_d(option = "magma") +
    facet_wrap(~temperature, scales = "free") +
    theme_bw() +
    scale_x_continuous(breaks = seq(min(dd$day), max(dd$day), by = 1)) +
    ggtitle(paste("Species:", species_name))
  
  # Save the ggplot object in the list§
  ggplot_objects[[species_name]] <- ggplot_object
}

# Access a specific plot, for example, Dexiostoma
dexiostoma_plot <- ggplot_objects[["Dexiostoma"]]
colpidium_plot <- ggplot_objects[["Colpidium"]]
loxocephalus_plot <- ggplot_objects[["Loxocephalus"]]
paramecium_plot <- ggplot_objects[["Paramecium"]]
spirostomum_plot <- ggplot_objects[["Spirostomum"]]
tetrahymena_plot <- ggplot_objects[["Tetrahymena"]]
# Display the plot
dexiostoma_plot # probably we can go to day 6 for Dexio
colpidium_plot
loxocephalus_plot
paramecium_plot
spirostomum_plot
tetrahymena_plot

# Create ggplot object
plot_regression_day <- ggplot(dd, aes(x = day, y = log(density + 1), color = as.factor(nutrients))) +
  geom_point(data = dd, aes(x = day, y = log(density + 1), color = as.factor(nutrients))) +
  geom_smooth(data = dd, method = "loess", se = FALSE) +
  geom_smooth(data = dd %>% filter(day <= 5), method = "lm", se = FALSE, linetype = "dashed", color = "black", aes(group = nutrients, color = nutrients)) +  # Add LM line
  scale_color_viridis_d(option = "magma") +
  facet_grid(temperature~ species, scales = "free") +
  theme_bw() +
  scale_x_continuous(breaks = seq(min(dd$day), max(dd$day), by = 1)) 




# Fit GAMS to calculated growth rate to estimate response surface

new_data <- expand_grid(temperature = seq(18, 28, by= 0.02),
                        nutrients = seq(1, 5, by= 0.01))


nested_gams <- dd_r_K %>% 
  nest(cols =-species) %>% 
  mutate(
    gams = map(cols, ~ gam(GrowthRate ~ ti(temperature) + ti(nutrients) + te(temperature, nutrients, k = c(5, 5)),
                           data = .x,
                           method = "REML")),
    predicted = map(gams, ~ predict(.x, newdata = new_data))
  )

# Get gams coefficients
coeff <- nested_gams %>% 
  mutate(coefs = map(gams, tidy, conf.int = TRUE)) %>% 
  unnest(coefs) %>% 
  dplyr::select(c('species', 'term', 'p.value'))



# Get gams glance
nested_gams %>% 
  mutate(results = map(gams, glance), 
         R.square = map_dbl(gams, ~ summary(.)$r.sq))  %>% 
  unnest(results) 


# Creating the dataset with the 2 columns as described above
predicted <- nested_gams %>% unnest(predicted)
rates <- cbind(new_data, predicted[, c(1,4)])
rates <- rates %>%
  relocate(species, temperature, nutrients, predicted)



# Create a list to store ggplot objects
surface_plots <- list()

# Loop through each species and draw GAM surface
for (i in seq_along(nested_gams$gams)) {
  species_name <- nested_gams$species[i]
  gam_model <- nested_gams$gams[[i]]
  
  # Create a new plot for each species
  plot <- ggplot(rates %>% filter(species == species_name),
                 aes(x = temperature, y = nutrients, z = predicted)) +
    geom_tile(aes(fill = predicted), color = "white") +
    geom_contour(color = "black", linetype = "solid", breaks = seq(-1, 1, by = 0.2)) +
    scale_fill_viridis_c() +
    labs(title = paste("GAM Surface for", species_name),
         fill = "Predicted Growth Rate") +
    theme_bw() +
    scale_x_continuous(breaks = c(18, 21, 24, 26, 28))
  
  # Save the plot in the list
  surface_plots[[species_name]] <- plot
  
}

# Access a specific plot, for example, Dexiostoma
surface_dexiostoma <- surface_plots[["Dexiostoma"]]
surface_colpidium <- surface_plots[["Colpidium"]]
surface_loxocephalus <- surface_plots[["Loxocephalus"]]
surface_paramecium <- surface_plots[["Paramecium"]]
surface_spirostotum <- surface_plots[["Spirostomum"]]
surface_tetrahymena <- surface_plots[["Tetrahymena"]]


# Checking predictions
dexiostoma_plot + surface_dexiostoma
colpidium_plot + surface_colpidium
loxocephalus_plot + surface_loxocephalus
paramecium_plot + surface_paramecium ### Paramecium seems to take longer to get to stable phase. Maybe do regression day 8
spirostomum_plot + surface_spirostotum ### Remove day 11 nutrients 2 temperature 28 
tetrahymena_plot + surface_tetrahymena


### Calculation Potential Response Diversity

refs2 <- crossing(temperature = seq(from = 18, to = 28, length.out = 10),
                  nutrients = seq(from = 1, to = 5, length.out = 10))


# Create all possible communities with 2 species

all_species <- c("Colpidium", "Dexiostoma", "Loxocephalus", "Paramecium", "Spirostomum", "Tetrahymena")
# Create all possible combinations of 2 species names
all_combinations_2species <- combn(c("Colpidium", "Dexiostoma", "Loxocephalus", "Paramecium", "Spirostomum", "Tetrahymena"), 2, simplify = TRUE)



# Calculate potential response diversity for all communities richness = 2


# Use lapply to iterate through each pair
# get_potential_RD <- function(species_names, nested_gams, E1_var, E2_var)  
# specify variables names for the "get_potential_RD" function
E1_var <- "temperature"
E2_var <- "nutrients"
result_list_2species <- lapply(1:ncol(all_combinations_2species), function(i) {
  E1_var <- "temperature"
  E2_var <- "nutrients"
  species_names <- all_combinations_2species[, i]
  result <- get_potential_RD(species_names, nested_gams, E1_var, E2_var)
  return(result)
})

# Combine all data frames into a single data frame
final_result_2species <- do.call(rbind, result_list_2species)

# Add a new column "composition" with the initial letters of species1 and species2
final_result_2species <- final_result_2species %>%
  mutate(composition = map_chr(row_number(), ~ paste0(substr(all_combinations_2species[, .], 1, 1), collapse = "")))%>% 
  mutate(richness = 2)
# Print the updated final result
print(final_result_2species)




### Calculating Response Diversity for all community compositions with richness = 3


# Generate all possible combinations of 3 species names
combinations_3species <- combn(all_species, 3, simplify = TRUE)

# specify variables names for the "get_potential_RD" function
E1_var <- "temperature"
E2_var <- "nutrients"
# Use lapply to iterate through each pair
result_list_3species <- lapply(1:ncol(combinations_3species), function(i) {
  species_names <- combinations_3species[, i]
  result <- get_potential_RD(species_names, nested_gams, E1_var, E2_var)
  return(result)
})


# Combine all data frames into a single data frame
result_list_3species <- do.call(rbind, result_list_3species)

# Add a new column "composition" with the initial letters of species1 and species2
result_list_3species <- result_list_3species %>%
  mutate(composition = map_chr(strsplit(species, "_"), ~ paste0(substr(., 1, 1), collapse = "")))%>% 
  mutate(richness = 3)
# Print the updated final result
print(result_list_3species)




# Calculating Response Diversity for all community compositions with richness = 4

# Generate all possible combinations of 3 species names
combinations_4species <- combn(all_species, 4, simplify = TRUE)

# specify variables names for the "get_potential_RD" function
E1_var <- "temperature"
E2_var <- "nutrients"
# Use lapply to iterate through each pair
result_list_4species <- lapply(1:ncol(combinations_4species), function(i) {
  species_names <- combinations_4species[, i]
  result <- get_potential_RD(species_names, nested_gams, E1_var, E2_var)
  return(result)
})


# Combine all data frames into a single data frame
result_list_4species <- do.call(rbind, result_list_4species)

# Add a new column "composition" with the initial letters of species1 and species2
result_list_4species <- result_list_4species %>%
  mutate(composition = map_chr(strsplit(species, "_"), ~ paste0(substr(., 1, 1), collapse = ""))) %>% 
  mutate(richness = 4)
# Print the updated final result
print(result_list_4species)


# Merging data frames and getting community compositions with the highest, lowest and intermediates potential response diversity 
# for each diversity level. 

# merge data frames
potential_RD <- rbind(final_result_2species, result_list_3species, result_list_4species)

potential_divergence_plot <- potential_RD %>% 
  ggplot(aes(x = richness, y = divergence)) +
  geom_point(size = 1.5) +
  geom_text(aes(label = composition), vjust = -0.5, size = 4) +
  theme_bw() + scale_x_continuous(breaks = c(2, 3, 4))  +
  labs(tag = "(b)")

potential_dissimilarity_plot <- potential_RD %>% 
  ggplot(aes(x = richness, y = dissimilarity)) +
  geom_point(size = 1.5) +
  geom_text(aes(label = composition), vjust = -0.5, size = 4) +
  theme_bw() + scale_x_continuous(breaks = c(2, 3, 4)) +
  labs(tag = "(a)")

potential_dissimilarity_plot + potential_divergence_plot



### data wrangling, probably not useful.
reshaped_potential_RD <- potential_RD %>%
  pivot_longer(cols = c(dissimilarity, divergence),
               names_to = "response_diversity",
               values_to = "values")



# Find rows with the highest and lowest values of response diversity
highest_values <- reshaped_potential_RD %>%
  group_by(richness, response_diversity) %>%
  filter(values == max(values)) %>%
  ungroup() %>% 
  mutate(potential_RD = "high")


lowest_values <- reshaped_potential_RD %>%
  group_by(richness, response_diversity) %>%
  filter(values == min(values)) %>%
  ungroup() %>% 
  mutate(potential_RD = "low")


mean_values <- reshaped_potential_RD %>%
  group_by(richness, response_diversity) %>%
  summarize(mean_value = mean(values))

# Find the rows with response diversity values closest to the median
mean_compositions <- reshaped_potential_RD %>%
  group_by(richness, response_diversity) %>%
  mutate(diff_from_mean = abs(values - mean_values$mean_value)) %>%
  filter(diff_from_mean == min(diff_from_mean)) %>%
  ungroup() %>% 
  select(species, composition, richness, response_diversity, values) %>% 
  mutate(potential_RD = "medium")


ranked_potential_RD_df <- rbind(highest_values, median_compositions, lowest_values)



#####################################################
####### REMEBER. #######################
## This changed compared to the .rmd version. This is the longer version analysing every T fluctuation in isolation



### Create a data set with environmental conditions that we may use in the experiment to see if RD calculated from known direction of environmental change 
### is different from potential RD, or if we get a good estimate.
### 3 fluctuating temperature x 3 fixed nutrients = 9 treatments

### T1 = fluctuating between 18 and 22 every 72h by 4 degree Celsius
# Create a sequence of time points (one measurement per day)
num_days = 60
time_points <- seq(0, num_days, by = 1)

# Define the combinations of nutrients for each time series
combinations <- list(
  combination_1 = 1,
  combination_2 = 3,
  combination_3 = 5
)

# Initialize an empty list to store the data frames for each combination
env_data_list <- list()

# Loop through each combination
for (i in seq_along(combinations)) {
  # Simulate temperature (fluctuating by +/- 2 degrees every 3 days)
  temperature <- 20 + 2 * sin(2 * pi * time_points / 3)  # Fluctuation every 3 days
  temperature <- pmin(pmax(temperature, 18), 22)  # Ensure temperature stays between 18 and 22
  
  # Create a data frame for the current combination
  env_data <- data.frame(
    time = time_points,  # Time in days
    temperature = temperature,
    nutrients = combinations[[i]],  # Set nutrients for the current combination
    combination = i  # Add combination identifier
  )
  
  # Append the data frame to the list
  env_data_list[[i]] <- env_data
}

# Combine the data frames for each combination into a single data frame
env_data_1 <- do.call(rbind, env_data_list)



# Plot temperature and nutrient level over time
# Plot temperature over time
p1 <- ggplot(env_data_1, aes(x = time, y = temperature)) +
  geom_line(color = "red") +
  labs(x = "Time (days)", y = "Temperature (°C)") +
  ggtitle("Temperature Over Time") +
  facet_grid(~combination) +
  ylim(18, 28) +
  theme_bw()

# Plot nutrient level over time
p2 <- ggplot(env_data_1, aes(x = time, y = nutrients)) +
  geom_line(color = "blue") +
  labs(x = "Time (days)", y = "Nutrient Level") +
  ggtitle("Nutrient Level Over Time") +
  facet_grid(~combination) +
  ylim(1, 5) +
  theme_bw()


env_data_1_plot <- p1 / p2



### T2 = fluctuating between 22 and 26 every 72h by 4 degree Celsius
num_days = 60
time_points <- seq(0, num_days, by = 1)

# Define the combinations of nutrients for each time series
combinations <- list(
  combination_1 = 1,
  combination_2 = 3,
  combination_3 = 5
)

# Initialize an empty list to store the data frames for each combination
env_data_list <- list()

# Loop through each combination
for (i in seq_along(combinations)) {
  # Simulate temperature (fluctuating by +/- 2 degrees every 3 days)
  temperature <- 24 + 2 * sin(2 * pi * time_points / 3)  # Fluctuation every 3 days
  temperature <- pmin(pmax(temperature, 22), 26)  # Ensure temperature stays between 18 and 22
  
  # Create a data frame for the current combination
  env_data <- data.frame(
    time = time_points,  # Time in days
    temperature = temperature,
    nutrients = combinations[[i]],  # Set nutrients for the current combination
    combination = i  # Add combination identifier
  )
  
  # Append the data frame to the list
  env_data_list[[i]] <- env_data
}

# Combine the data frames for each combination into a single data frame
env_data_2 <- do.call(rbind, env_data_list)

# Plot temperature and nutrient level over time
# Plot temperature over time
p1 <- ggplot(env_data_2, aes(x = time, y = temperature)) +
  geom_line(color = "red") +
  labs(x = "Time (days)", y = "Temperature (°C)") +
  ggtitle("Temperature Over Time") +
  facet_grid(~combination) +
  ylim(18, 28) +
  theme_bw()

# Plot nutrient level over time
p2 <- ggplot(env_data_2, aes(x = time, y = nutrients)) +
  geom_line(color = "blue") +
  labs(x = "Time (days)", y = "Nutrient Level") +
  ggtitle("Nutrient Level Over Time") +
  facet_grid(~combination) +
  ylim(1, 5) +
  theme_bw()


env_data_2_plot <- p1 / p2



### T3 = fluctuating between 24 and 28 every 72h by 4 degree Celsius
num_days = 60
time_points <- seq(0, num_days, by = 1)

# Define the combinations of nutrients for each time series
combinations <- list(
  combination_1 = 1,
  combination_2 = 3,
  combination_3 = 5
)

# Initialize an empty list to store the data frames for each combination
env_data_list <- list()

# Loop through each combination
for (i in seq_along(combinations)) {
  # Simulate temperature (fluctuating by +/- 2 degrees every 3 days)
  temperature <- 26 + 2 * sin(2 * pi * time_points / 3)  # Fluctuation every 3 days
  temperature <- pmin(pmax(temperature, 24), 28)  # Ensure temperature stays between 18 and 22
  
  # Create a data frame for the current combination
  env_data <- data.frame(
    time = time_points,  # Time in days
    temperature = temperature,
    nutrients = combinations[[i]],  # Set nutrients for the current combination
    combination = i  # Add combination identifier
  )
  
  # Append the data frame to the list
  env_data_list[[i]] <- env_data
}

# Combine the data frames for each combination into a single data frame
env_data_3 <- do.call(rbind, env_data_list)

# Plot temperature and nutrient level over time
# Plot temperature over time
p1 <- ggplot(env_data_3, aes(x = time, y = temperature)) +
  geom_line(color = "red") +
  labs(x = "Time (days)", y = "Temperature (°C)") +
  ggtitle("Temperature Over Time") +
  facet_grid(~combination) +
  ylim(18, 28) +
  theme_bw()

# Plot nutrient level over time
p2 <- ggplot(env_data_3, aes(x = time, y = nutrients)) +
  geom_line(color = "blue") +
  labs(x = "Time (days)", y = "Nutrient Level") +
  ggtitle("Nutrient Level Over Time") +
  facet_grid(~combination) +
  ylim(1, 5) +
  theme_bw()


env_data_3_plot <- p1 / p2


### Try to merge the 3 env_data
env_data_1 <- env_data_1 %>% mutate(temperature_treatment = "18_22")
env_data_2 <- env_data_2 %>% mutate(temperature_treatment = "22_26")
env_data_3 <- env_data_3 %>% mutate(temperature_treatment = "24_28")
env_data <- rbind(env_data_1, env_data_2, env_data_3)


E1_var <- "temperature"
E2_var <- "nutrients"
env_data <- env_data
result_list_2species_RD <- lapply(1:ncol(all_combinations_2species), function(i) {
  E1_var <- "temperature"
  E2_var <- "nutrients"
  species_names <- all_combinations_2species[, i]
  result <- get_RD(species_names, nested_gams, E1_var, E2_var)
  return(result)
})


# Combine all data frames into a single data frame
RD_2species <- do.call(rbind, result_list_2species_RD)

# Add a new column "composition" with the initial letters of species1 and species2
RD_2species <- RD_2species %>%
  mutate(composition = str_split(species, "_") %>% 
           map_chr(~ paste0(substr(.x, 1, 1), collapse = ""))) %>%   mutate(richness = 2)
# Print the updated final result
print(RD_2species)


RD_2species <- RD_2species %>% rename(temperature = E1_ref,
                         nutrients = E2_ref) %>% 
  arrange(nutrients) %>% 
  mutate(combination = case_when(
    nutrients == 1 ~ 1,
    nutrients == 3 ~ 2,
    nutrients == 5 ~ 3
  ))


# Plot boxplots for divergence by composition using different colors from the Viridis palette
dissimilarity_2sp_plot <- ggplot(RD_2species, aes(x = composition, y = dissimilarity, fill = composition)) +
  geom_boxplot(alpha = 0.7) +  # Boxplots
  geom_jitter(width = 0.2, alpha = 0.7) +  # Individual points with jitter
  labs(x = "Composition", y = "Dissimilarity") +
  #facet_grid(~combination) +
  scale_fill_viridis_d() +  # Use the Viridis palette
  theme_bw()


divergence_2sp_plot <-  ggplot(RD_2species_1, aes(x = composition, y = divergence, fill = composition)) +
  geom_boxplot(alpha = 0.7) +  # Boxplots
  geom_jitter(width = 0.2, alpha = 0.7) +  # Individual points with jitter
  labs(x = "Composition", y = "Divergence") +
  #facet_grid(~combination) +
  scale_fill_viridis_d() +  # Use the Viridis palette
  theme_bw()

RD_T18_22_2sp  <- dissimilarity_2sp_plot + divergence_2sp_plot

### Calculate the composition with the highest response divergence across the 3 combinations
summary_RD_2spp_1 <- RD_2species_1 %>% dplyr::group_by(composition) %>% 
  mutate(mean_div = mean(divergence)) %>% 
  select(composition, mean_div, richness) %>% 
  distinct(composition, .keep_all = TRUE)

summary_RD_2spp_1 %>% 
  ggplot(aes(x = richness, y = mean_div)) +
  geom_point(size = 1.5) +
  geom_text(aes(label = composition), vjust = -0.5, size = 4) +
  theme_bw() + scale_x_continuous(breaks = c(2, 3, 4)) 





### T2 = fluctuating between 22 and 26 every 72h by 4 degree Celsius
# Create a sequence of time points (one measurement per day)
time_points <- seq(0, num_days, by = 1)

# Define the combinations of nutrients for each time series
combinations <- list(
  combination_1 = 1,
  combination_2 = 3,
  combination_3 = 5
)

# Initialize an empty list to store the data frames for each combination
env_data_list <- list()

# Loop through each combination
for (i in seq_along(combinations)) {
  # Simulate temperature (fluctuating by +/- 2 degrees every 3 days)
  temperature <- 24 + 2 * sin(2 * pi * time_points / 3)  # Fluctuation every 3 days
  temperature <- pmin(pmax(temperature, 22), 26)  # Ensure temperature stays between 18 and 22
  
  # Create a data frame for the current combination
  env_data <- data.frame(
    time = time_points,  # Time in days
    temperature = temperature,
    nutrients = combinations[[i]],  # Set nutrients for the current combination
    combination = i  # Add combination identifier
  )
  
  # Append the data frame to the list
  env_data_list[[i]] <- env_data
}

# Combine the data frames for each combination into a single data frame
env_data_2 <- do.call(rbind, env_data_list)



# Plot temperature and nutrient level over time
# Plot temperature over time
p1 <- ggplot(env_data_2, aes(x = time, y = temperature)) +
  geom_line(color = "red") +
  labs(x = "Time (days)", y = "Temperature (°C)") +
  ggtitle("Temperature Over Time") +
  facet_grid(~combination) +
  ylim(18, 28) +
  theme_bw()

# Plot nutrient level over time
p2 <- ggplot(env_data_2, aes(x = time, y = nutrients)) +
  geom_line(color = "blue") +
  labs(x = "Time (days)", y = "Nutrient Level") +
  ggtitle("Nutrient Level Over Time") +
  facet_grid(~combination) +
  ylim(1, 5) +
  theme_bw()


p1 / p2


E1_var <- "temperature"
E2_var <- "nutrients"
env_data <- env_data_2
result_list_2species_RD_2 <- lapply(1:ncol(all_combinations_2species), function(i) {
  E1_var <- "temperature"
  E2_var <- "nutrients"
  species_names <- all_combinations_2species[, i]
  result <- get_RD(species_names, nested_gams, E1_var, E2_var)
  return(result)
})


# Combine all data frames into a single data frame
RD_2species_2 <- do.call(rbind, result_list_2species_RD_2)

# Add a new column "composition" with the initial letters of species1 and species2
RD_2species_2 <- RD_2species_2 %>%
  mutate(composition = str_split(species, "_") %>% 
           map_chr(~ paste0(substr(.x, 1, 1), collapse = ""))) %>%   mutate(richness = 2)
# Print the updated final result
print(RD_2species_2)


RD_2species_2 <- RD_2species_2 %>% rename(temperature = E1_ref,
                                          nutrients = E2_ref) %>% 
  arrange(nutrients) %>% 
  mutate(combination = case_when(
    nutrients == 1 ~ 1,
    nutrients == 3 ~ 2,
    nutrients == 5 ~ 3
  ))


# Plot boxplots for divergence by composition using different colors from the Viridis palette
dissimilarity_2sp_plot_2 <- ggplot(RD_2species_2, aes(x = composition, y = dissimilarity, fill = composition)) +
  geom_boxplot(alpha = 0.7) +  # Boxplots
  geom_jitter(width = 0.2, alpha = 0.7) +  # Individual points with jitter
  labs(x = "Composition", y = "Dissimilarity") +
  facet_grid(~combination) +
  scale_fill_viridis_d() +  # Use the Viridis palette
  theme_bw()


divergence_2sp_plot_2 <-  ggplot(RD_2species_2, aes(x = composition, y = divergence, fill = composition)) +
  geom_boxplot(alpha = 0.7) +  # Boxplots
  geom_jitter(width = 0.2, alpha = 0.7) +  # Individual points with jitter
  labs(x = "Composition", y = "Divergence") +
  facet_grid(~combination) +
  scale_fill_viridis_d() +  # Use the Viridis palette
  theme_bw()

RD_T22_26_2sp  <- dissimilarity_2sp_plot_2 + divergence_2sp_plot_2

### Calculate the composition with the highest response divergence across the 3 combinations
summary_RD_2spp_2 <- RD_2species_2 %>% dplyr::group_by(composition) %>% 
  mutate(mean_div = mean(divergence)) %>% 
  select(composition, mean_div, richness) %>% 
  distinct(composition, .keep_all = TRUE)

summary_RD_2spp_2 %>% 
  ggplot(aes(x = richness, y = mean_div)) +
  geom_point(size = 1.5) +
  geom_text(aes(label = composition), vjust = -0.5, size = 4) +
  theme_bw() + scale_x_continuous(breaks = c(2, 3, 4)) 




### T3 = fluctuating between 24 and 28 every 72h by 4 degree Celsius
# Create a sequence of time points (one measurement per day)
time_points <- seq(0, num_days, by = 1)

# Define the combinations of nutrients for each time series
combinations <- list(
  combination_1 = 1,
  combination_2 = 3,
  combination_3 = 5
)

# Initialize an empty list to store the data frames for each combination
env_data_list <- list()

# Loop through each combination
for (i in seq_along(combinations)) {
  # Simulate temperature (fluctuating by +/- 2 degrees every 3 days)
  temperature <- 26 + 2 * sin(2 * pi * time_points / 3)  # Fluctuation every 3 days
  temperature <- pmin(pmax(temperature, 24), 28)  # Ensure temperature stays between 18 and 22
  
  # Create a data frame for the current combination
  env_data <- data.frame(
    time = time_points,  # Time in days
    temperature = temperature,
    nutrients = combinations[[i]],  # Set nutrients for the current combination
    combination = i  # Add combination identifier
  )
  
  # Append the data frame to the list
  env_data_list[[i]] <- env_data
}

# Combine the data frames for each combination into a single data frame
env_data_3 <- do.call(rbind, env_data_list)



# Plot temperature and nutrient level over time
# Plot temperature over time
p1 <- ggplot(env_data_3, aes(x = time, y = temperature)) +
  geom_line(color = "red") +
  labs(x = "Time (days)", y = "Temperature (°C)") +
  ggtitle("Temperature Over Time") +
  facet_grid(~combination) +
  ylim(18, 28) +
  theme_bw()

# Plot nutrient level over time
p2 <- ggplot(env_data_3, aes(x = time, y = nutrients)) +
  geom_line(color = "blue") +
  labs(x = "Time (days)", y = "Nutrient Level") +
  ggtitle("Nutrient Level Over Time") +
  facet_grid(~combination) +
  ylim(1, 5) +
  theme_bw()


p1 / p2


E1_var <- "temperature"
E2_var <- "nutrients"
env_data <- env_data_3
result_list_2species_RD_3 <- lapply(1:ncol(all_combinations_2species), function(i) {
  E1_var <- "temperature"
  E2_var <- "nutrients"
  species_names <- all_combinations_2species[, i]
  result <- get_RD(species_names, nested_gams, E1_var, E2_var)
  return(result)
})


# Combine all data frames into a single data frame
RD_2species_3 <- do.call(rbind, result_list_2species_RD_3)

# Add a new column "composition" with the initial letters of species1 and species2
RD_2species_3 <- RD_2species_3 %>%
  mutate(composition = str_split(species, "_") %>% 
           map_chr(~ paste0(substr(.x, 1, 1), collapse = ""))) %>%   mutate(richness = 2)
# Print the updated final result
print(RD_2species_3)


RD_2species_3 <- RD_2species_3 %>% rename(temperature = E1_ref,
                                          nutrients = E2_ref) %>% 
  arrange(nutrients) %>% 
  mutate(combination = case_when(
    nutrients == 1 ~ 1,
    nutrients == 3 ~ 2,
    nutrients == 5 ~ 3
  ))


# Plot boxplots for divergence by composition using different colors from the Viridis palette
dissimilarity_2sp_plot_3 <- ggplot(RD_2species_3, aes(x = composition, y = dissimilarity, fill = composition)) +
  geom_boxplot(alpha = 0.7) +  # Boxplots
  geom_jitter(width = 0.2, alpha = 0.7) +  # Individual points with jitter
  labs(x = "Composition", y = "Dissimilarity") +
  facet_grid(~combination) +
  scale_fill_viridis_d() +  # Use the Viridis palette
  theme_bw()


divergence_2sp_plot_3 <-  ggplot(RD_2species_3, aes(x = composition, y = divergence, fill = composition)) +
  geom_boxplot(alpha = 0.7) +  # Boxplots
  geom_jitter(width = 0.2, alpha = 0.7) +  # Individual points with jitter
  labs(x = "Composition", y = "Divergence") +
  facet_grid(~combination) +
  scale_fill_viridis_d() +  # Use the Viridis palette
  theme_bw()

RD_T24_28_2sp <- dissimilarity_2sp_plot_3 + divergence_2sp_plot_3



### Calculate the composition with the highest response divergence across the 3 combinations
summary_RD_2spp_3 <- RD_2species_3 %>% dplyr::group_by(composition) %>% 
  mutate(mean_div = mean(divergence)) %>% 
  select(composition, mean_div, richness) %>% 
  distinct(composition, .keep_all = TRUE)

summary_RD_2spp_3 %>% 
  ggplot(aes(x = richness, y = mean_div)) +
  geom_point(size = 1.5) +
  geom_text(aes(label = composition), vjust = -0.5, size = 4) +
  theme_bw() + scale_x_continuous(breaks = c(2, 3, 4)) 


summary_RD_2 <- rbind(summary_RD_2spp_1, summary_RD_2spp_2, summary_RD_2spp_3)
summary_RD_2 <- summary_RD_2 %>% 
  dplyr::group_by(composition) %>% 
  mutate(mean_div2 = mean(mean_div)) %>% 
  select(composition, mean_div2, richness) %>% 
  distinct(composition, .keep_all = TRUE)
summary_RD_2 %>% 
  ggplot(aes(x = richness, y = mean_div2)) +
  geom_point(size = 1.5) +
  geom_text(aes(label = composition), vjust = -0.5, size = 4) +
  theme_bw() + scale_x_continuous(breaks = c(2, 3, 4)) 



### Now same but with 3 species
#T1
E1_var <- "temperature"
E2_var <- "nutrients"
env_data <- env_data_1
result_list_3species_RD <- lapply(1:ncol(combinations_3species), function(i) {
  E1_var <- "temperature"
  E2_var <- "nutrients"
  species_names <- combinations_3species[, i]
  result <- get_RD(species_names, nested_gams, E1_var, E2_var)
  return(result)
})


# Combine all data frames into a single data frame
RD_3species_1 <- do.call(rbind, result_list_3species_RD)

# Add a new column "composition" with the initial letters of species1 and species2
RD_3species_1 <- RD_3species_1 %>%
  mutate(composition = str_split(species, "_") %>% 
           map_chr(~ paste0(substr(.x, 1, 1), collapse = ""))) %>%   mutate(richness = 3)
# Print the updated final result
print(RD_3species_1)


RD_3species_1 <- RD_3species_1 %>% rename(temperature = E1_ref,
                                          nutrients = E2_ref) %>% 
  arrange(nutrients) %>% 
  mutate(combination = case_when(
    nutrients == 1 ~ 1,
    nutrients == 3 ~ 2,
    nutrients == 5 ~ 3
  ))


# Plot boxplots for divergence by composition using different colors from the Viridis palette
dissimilarity_3sp_plot <- ggplot(RD_3species_1, aes(x = composition, y = dissimilarity, fill = composition)) +
  geom_boxplot(alpha = 0.7) +  # Boxplots
  geom_jitter(width = 0.2, alpha = 0.7) +  # Individual points with jitter
  labs(x = "Composition", y = "Dissimilarity") +
  facet_grid(~combination) +
  scale_fill_viridis_d() +  # Use the Viridis palette
  theme_bw()


divergence_3sp_plot <-  ggplot(RD_3species_1, aes(x = composition, y = divergence, fill = composition)) +
  geom_boxplot(alpha = 0.7) +  # Boxplots
  geom_jitter(width = 0.2, alpha = 0.7) +  # Individual points with jitter
  labs(x = "Composition", y = "Divergence") +
  facet_grid(~combination) +
  scale_fill_viridis_d() +  # Use the Viridis palette
  theme_bw()

RD_T18_22_3sp  <- dissimilarity_3sp_plot + divergence_3sp_plot

### Calculate the composition with the highest response divergence across the 3 combinations
summary_RD_3spp_1 <- RD_3species_1 %>% dplyr::group_by(composition) %>% 
  mutate(mean_div = mean(divergence)) %>% 
  select(composition, mean_div, richness) %>% 
  distinct(composition, .keep_all = TRUE)

summary_RD_3spp_1 %>% 
  ggplot(aes(x = richness, y = mean_div)) +
  geom_point(size = 1.5) +
  geom_text(aes(label = composition), vjust = -0.5, size = 4) +
  theme_bw() + scale_x_continuous(breaks = c(2, 3, 4)) 


#T2
E1_var <- "temperature"
E2_var <- "nutrients"
env_data <- env_data_2
result_list_3species_RD <- lapply(1:ncol(combinations_3species), function(i) {
  E1_var <- "temperature"
  E2_var <- "nutrients"
  species_names <- combinations_3species[, i]
  result <- get_RD(species_names, nested_gams, E1_var, E2_var)
  return(result)
})


# Combine all data frames into a single data frame
RD_3species_2 <- do.call(rbind, result_list_3species_RD)

# Add a new column "composition" with the initial letters of species1 and species2
RD_3species_2 <- RD_3species_2 %>%
  mutate(composition = str_split(species, "_") %>% 
           map_chr(~ paste0(substr(.x, 1, 1), collapse = ""))) %>%   mutate(richness = 3)
# Print the updated final result
print(RD_3species_2)


RD_3species_2 <- RD_3species_2 %>% rename(temperature = E1_ref,
                                          nutrients = E2_ref) %>% 
  arrange(nutrients) %>% 
  mutate(combination = case_when(
    nutrients == 1 ~ 1,
    nutrients == 3 ~ 2,
    nutrients == 5 ~ 3
  ))


# Plot boxplots for divergence by composition using different colors from the Viridis palette
dissimilarity_3sp_plot_2 <- ggplot(RD_3species_2, aes(x = composition, y = dissimilarity, fill = composition)) +
  geom_boxplot(alpha = 0.7) +  # Boxplots
  geom_jitter(width = 0.2, alpha = 0.7) +  # Individual points with jitter
  labs(x = "Composition", y = "Dissimilarity") +
  facet_grid(~combination) +
  scale_fill_viridis_d() +  # Use the Viridis palette
  theme_bw()


divergence_3sp_plot_2 <-  ggplot(RD_3species_2, aes(x = composition, y = divergence, fill = composition)) +
  geom_boxplot(alpha = 0.7) +  # Boxplots
  geom_jitter(width = 0.2, alpha = 0.7) +  # Individual points with jitter
  labs(x = "Composition", y = "Divergence") +
  facet_grid(~combination) +
  scale_fill_viridis_d() +  # Use the Viridis palette
  theme_bw()

RD_T22_26_3sp  <- dissimilarity_3sp_plot_2 + divergence_3sp_plot_2


### Calculate the composition with the highest response divergence across the 3 combinations
summary_RD_3spp_2 <- RD_3species_2 %>% dplyr::group_by(composition) %>% 
  mutate(mean_div = mean(divergence)) %>% 
  select(composition, mean_div, richness) %>% 
  distinct(composition, .keep_all = TRUE)

summary_RD_3spp_2 %>% 
  ggplot(aes(x = richness, y = mean_div)) +
  geom_point(size = 1.5) +
  geom_text(aes(label = composition), vjust = -0.5, size = 4) +
  theme_bw() + scale_x_continuous(breaks = c(2, 3, 4)) 


#T3
E1_var <- "temperature"
E2_var <- "nutrients"
env_data <- env_data_3
result_list_3species_RD <- lapply(1:ncol(combinations_3species), function(i) {
  E1_var <- "temperature"
  E2_var <- "nutrients"
  species_names <- combinations_3species[, i]
  result <- get_RD(species_names, nested_gams, E1_var, E2_var)
  return(result)
})


# Combine all data frames into a single data frame
RD_3species_3 <- do.call(rbind, result_list_3species_RD)

# Add a new column "composition" with the initial letters of species1 and species2
RD_3species_3 <- RD_3species_3 %>%
  mutate(composition = str_split(species, "_") %>% 
           map_chr(~ paste0(substr(.x, 1, 1), collapse = ""))) %>%   mutate(richness = 3)
# Print the updated final result
print(RD_3species_3)


RD_3species_3 <- RD_3species_3 %>% rename(temperature = E1_ref,
                                          nutrients = E2_ref) %>% 
  arrange(nutrients) %>% 
  mutate(combination = case_when(
    nutrients == 1 ~ 1,
    nutrients == 3 ~ 2,
    nutrients == 5 ~ 3
  ))


# Plot boxplots for divergence by composition using different colors from the Viridis palette
dissimilarity_3sp_plot_3 <- ggplot(RD_3species_3, aes(x = composition, y = dissimilarity, fill = composition)) +
  geom_boxplot(alpha = 0.7) +  # Boxplots
  geom_jitter(width = 0.2, alpha = 0.7) +  # Individual points with jitter
  labs(x = "Composition", y = "Dissimilarity") +
  facet_grid(~combination) +
  scale_fill_viridis_d() +  # Use the Viridis palette
  theme_bw()


divergence_3sp_plot_3 <-  ggplot(RD_3species_3, aes(x = composition, y = divergence, fill = composition)) +
  geom_boxplot(alpha = 0.7) +  # Boxplots
  geom_jitter(width = 0.2, alpha = 0.7) +  # Individual points with jitter
  labs(x = "Composition", y = "Divergence") +
  facet_grid(~combination) +
  scale_fill_viridis_d() +  # Use the Viridis palette
  theme_bw()

RD_T24_28_3sp  <- dissimilarity_3sp_plot_3 + divergence_3sp_plot_3


### Calculate the composition with the highest response divergence across the 3 combinations
summary_RD_3spp_3 <- RD_3species_3 %>% dplyr::group_by(composition) %>% 
  mutate(mean_div = mean(divergence)) %>% 
  select(composition, mean_div, richness) %>% 
  distinct(composition, .keep_all = TRUE)

summary_RD_3spp_3 %>% 
  ggplot(aes(x = richness, y = mean_div)) +
  geom_point(size = 1.5) +
  geom_text(aes(label = composition), vjust = -0.5, size = 4) +
  theme_bw() + scale_x_continuous(breaks = c(2, 3, 4)) 



summary_RD_3 <- rbind(summary_RD_3spp_1, summary_RD_3spp_2, summary_RD_3spp_3)
summary_RD_3 <- summary_RD_3 %>% 
  dplyr::group_by(composition) %>% 
  mutate(mean_div2 = mean(mean_div)) %>% 
  select(composition, mean_div2, richness) %>% 
  distinct(composition, .keep_all = TRUE)
summary_RD_3 %>% 
  ggplot(aes(x = richness, y = mean_div2)) +
  geom_point(size = 1.5) +
  geom_text(aes(label = composition), vjust = -0.5, size = 4) +
  theme_bw() + scale_x_continuous(breaks = c(2, 3, 4)) 


### Now same but with 4 species
#T1
E1_var <- "temperature"
E2_var <- "nutrients"
env_data <- env_data_1
result_list_4species_RD <- lapply(1:ncol(combinations_4species), function(i) {
  E1_var <- "temperature"
  E2_var <- "nutrients"
  species_names <- combinations_4species[, i]
  result <- get_RD(species_names, nested_gams, E1_var, E2_var)
  return(result)
})


# Combine all data frames into a single data frame
RD_4species_1 <- do.call(rbind, result_list_4species_RD)

# Add a new column "composition" with the initial letters of species1 and species2
RD_4species_1 <- RD_4species_1 %>%
  mutate(composition = str_split(species, "_") %>% 
           map_chr(~ paste0(substr(.x, 1, 1), collapse = ""))) %>%   mutate(richness = 4)
# Print the updated final result
print(RD_4species_1)


RD_4species_1 <- RD_4species_1 %>% rename(temperature = E1_ref,
                                          nutrients = E2_ref) %>% 
  arrange(nutrients) %>% 
  mutate(combination = case_when(
    nutrients == 1 ~ 1,
    nutrients == 3 ~ 2,
    nutrients == 5 ~ 3
  ))


# Plot boxplots for divergence by composition using different colors from the Viridis palette
dissimilarity_4sp_plot <- ggplot(RD_4species_1, aes(x = composition, y = dissimilarity, fill = composition)) +
  geom_boxplot(alpha = 0.7) +  # Boxplots
  geom_jitter(width = 0.2, alpha = 0.7) +  # Individual points with jitter
  labs(x = "Composition", y = "Dissimilarity") +
  facet_grid(~combination) +
  scale_fill_viridis_d() +  # Use the Viridis palette
  theme_bw()


divergence_4sp_plot <-  ggplot(RD_4species_1, aes(x = composition, y = divergence, fill = composition)) +
  geom_boxplot(alpha = 0.7) +  # Boxplots
  geom_jitter(width = 0.2, alpha = 0.7) +  # Individual points with jitter
  labs(x = "Composition", y = "Divergence") +
  facet_grid(~combination) +
  scale_fill_viridis_d() +  # Use the Viridis palette
  theme_bw()

RD_T18_22_4sp  <- dissimilarity_4sp_plot + divergence_4sp_plot

### Calculate the composition with the highest response divergence across the 3 combinations
summary_RD_4spp_1 <- RD_4species_1 %>% dplyr::group_by(composition) %>% 
  mutate(mean_div = mean(divergence)) %>% 
  select(composition, mean_div, richness) %>% 
  distinct(composition, .keep_all = TRUE)

summary_RD_4spp_1 %>% 
  ggplot(aes(x = richness, y = mean_div)) +
  geom_point(size = 1.5) +
  geom_text(aes(label = composition), vjust = -0.5, size = 4) +
  theme_bw() + scale_x_continuous(breaks = c(2, 3, 4)) 

#T2
E1_var <- "temperature"
E2_var <- "nutrients"
env_data <- env_data_2
result_list_4species_RD <- lapply(1:ncol(combinations_4species), function(i) {
  E1_var <- "temperature"
  E2_var <- "nutrients"
  species_names <- combinations_4species[, i]
  result <- get_RD(species_names, nested_gams, E1_var, E2_var)
  return(result)
})


# Combine all data frames into a single data frame
RD_4species_2 <- do.call(rbind, result_list_4species_RD)

# Add a new column "composition" with the initial letters of species1 and species2
RD_4species_2 <- RD_4species_2 %>%
  mutate(composition = str_split(species, "_") %>% 
           map_chr(~ paste0(substr(.x, 1, 1), collapse = ""))) %>%   mutate(richness = 3)
# Print the updated final result
print(RD_4species_2)


RD_4species_2 <- RD_4species_2 %>% rename(temperature = E1_ref,
                                          nutrients = E2_ref) %>% 
  arrange(nutrients) %>% 
  mutate(combination = case_when(
    nutrients == 1 ~ 1,
    nutrients == 3 ~ 2,
    nutrients == 5 ~ 3
  ))


# Plot boxplots for divergence by composition using different colors from the Viridis palette
dissimilarity_4sp_plot_2 <- ggplot(RD_4species_2, aes(x = composition, y = dissimilarity, fill = composition)) +
  geom_boxplot(alpha = 0.7) +  # Boxplots
  geom_jitter(width = 0.2, alpha = 0.7) +  # Individual points with jitter
  labs(x = "Composition", y = "Dissimilarity") +
  facet_grid(~combination) +
  scale_fill_viridis_d() +  # Use the Viridis palette
  theme_bw()


divergence_4sp_plot_2 <-  ggplot(RD_4species_2, aes(x = composition, y = divergence, fill = composition)) +
  geom_boxplot(alpha = 0.7) +  # Boxplots
  geom_jitter(width = 0.2, alpha = 0.7) +  # Individual points with jitter
  labs(x = "Composition", y = "Divergence") +
  facet_grid(~combination) +
  scale_fill_viridis_d() +  # Use the Viridis palette
  theme_bw()

RD_T22_26_4sp  <- dissimilarity_4sp_plot_2 + divergence_4sp_plot_2

### Calculate the composition with the highest response divergence across the 3 combinations
summary_RD_4spp_2 <- RD_4species_2 %>% dplyr::group_by(composition) %>% 
  mutate(mean_div = mean(divergence)) %>% 
  select(composition, mean_div, richness) %>% 
  distinct(composition, .keep_all = TRUE)

summary_RD_4spp_2 %>% 
  ggplot(aes(x = richness, y = mean_div)) +
  geom_point(size = 1.5) +
  geom_text(aes(label = composition), vjust = -0.5, size = 4) +
  theme_bw() + scale_x_continuous(breaks = c(2, 3, 4)) 

#T3
E1_var <- "temperature"
E2_var <- "nutrients"
env_data <- env_data_3
result_list_4species_RD <- lapply(1:ncol(combinations_4species), function(i) {
  E1_var <- "temperature"
  E2_var <- "nutrients"
  species_names <- combinations_4species[, i]
  result <- get_RD(species_names, nested_gams, E1_var, E2_var)
  return(result)
})


# Combine all data frames into a single data frame
RD_4species_3 <- do.call(rbind, result_list_4species_RD)

# Add a new column "composition" with the initial letters of species1 and species2
RD_4species_3 <- RD_4species_3 %>%
  mutate(composition = str_split(species, "_") %>% 
           map_chr(~ paste0(substr(.x, 1, 1), collapse = ""))) %>%   mutate(richness = 3)
# Print the updated final result
print(RD_4species_3)


RD_4species_3 <- RD_4species_3 %>% rename(temperature = E1_ref,
                                          nutrients = E2_ref) %>% 
  arrange(nutrients) %>% 
  mutate(combination = case_when(
    nutrients == 1 ~ 1,
    nutrients == 3 ~ 2,
    nutrients == 5 ~ 3
  ))


# Plot boxplots for divergence by composition using different colors from the Viridis palette
dissimilarity_4sp_plot_3 <- ggplot(RD_4species_3, aes(x = composition, y = dissimilarity, fill = composition)) +
  geom_boxplot(alpha = 0.7) +  # Boxplots
  geom_jitter(width = 0.2, alpha = 0.7) +  # Individual points with jitter
  labs(x = "Composition", y = "Dissimilarity") +
  facet_grid(~combination) +
  scale_fill_viridis_d() +  # Use the Viridis palette
  theme_bw()


divergence_4sp_plot_3 <-  ggplot(RD_4species_3, aes(x = composition, y = divergence, fill = composition)) +
  geom_boxplot(alpha = 0.7) +  # Boxplots
  geom_jitter(width = 0.2, alpha = 0.7) +  # Individual points with jitter
  labs(x = "Composition", y = "Divergence") +
  facet_grid(~combination) +
  scale_fill_viridis_d() +  # Use the Viridis palette
  theme_bw()

RD_T24_28_4sp  <- dissimilarity_4sp_plot_3 + divergence_4sp_plot_3

### Calculate the composition with the highest response divergence across the 3 combinations
summary_RD_4spp_3 <- RD_4species_3 %>% dplyr::group_by(composition) %>% 
  mutate(mean_div = mean(divergence)) %>% 
  select(composition, mean_div, richness) %>% 
  distinct(composition, .keep_all = TRUE)

summary_RD_4spp_3 %>% 
  ggplot(aes(x = richness, y = mean_div)) +
  geom_point(size = 1.5) +
  geom_text(aes(label = composition), vjust = -0.5, size = 4) +
  theme_bw() + scale_x_continuous(breaks = c(2, 3, 4)) 


summary_RD_4 <- rbind(summary_RD_4spp_1, summary_RD_4spp_2, summary_RD_4spp_3)
summary_RD_4 <- summary_RD_4 %>% 
  dplyr::group_by(composition) %>% 
  mutate(mean_div2 = mean(mean_div)) %>% 
  select(composition, mean_div2, richness) %>% 
  distinct(composition, .keep_all = TRUE)


summary_RD_4 %>% 
  ggplot(aes(x = richness, y = mean_div2)) +
  geom_point(size = 1.5) +
  geom_text(aes(label = composition), vjust = -0.5, size = 4) +
  theme_bw() + scale_x_continuous(breaks = c(2, 3, 4)) 


summary_RD <- rbind(summary_RD_2, summary_RD_3, summary_RD_4)

summary_RD %>% 
  ggplot(aes(x = richness, y = mean_div2)) +
  geom_point(size = 1.5) +
  geom_text(aes(label = composition), vjust = -0.5, size = 4) +
  theme_bw() + scale_x_continuous(breaks = c(2, 3, 4)) +
  labs(y = "Mean Divergence")
