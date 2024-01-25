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

### Dummy try to fig GAMs with tensor product to simulated data and calculate r and K
# Install and load required packages
# install.packages(c("mgcv", "ggplot2"))
library(mgcv)
library(ggplot2)

# Simulated data with 6 species, 5 temperature levels, and 5 nutrients levels
set.seed(123)

# Number of observations per combination
n <- 100

# Create a data frame with 6 species
species <- rep(letters[1:6], each = 5 * n)

# Create a data frame with 5 temperature levels and 5 nutrients levels
temp <- rep(1:5, each = n)
nutrients <- rep(1:5, times = n)

# Simulate response variable for each species
response <- rnorm(length(species), mean = 10 + 2 * temp + 3 * nutrients, sd = 5)

# Create a data frame with numeric variables
data <- data.frame(Species = factor(species),
                   Temp = rep(temp, times = 6),
                   Nutrients = rep(nutrients, times = 6),
                   Response = response)

# Fit GAM with tensor product for the factorial design
gam_model <- gam(Response ~ te(Temp, Nutrients, by = Species), data = data)

# Summary of the GAM model
summary(gam_model)

# Predict response surface for visualization
new_data <- expand.grid(Species = levels(data$Species),
                        Temp = seq(min(data$Temp), max(data$Temp), length.out = 100),
                        Nutrients = seq(min(data$Nutrients), max(data$Nutrients), length.out = 100))

predictions <- predict(gam_model, newdata = new_data, type = "response")

# Plot the response surface for one species (change 'a' to the desired species)
ggplot(data, aes(x = Temp, y = Nutrients, z = Response, color = Species)) +
  geom_tile(aes(fill = Response), color = "white") +
  geom_contour(data = as.data.frame(new_data), aes(x = Temp, y = Nutrients, z = predictions), color = "black") +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "Response Surface for All Species") +
  theme_minimal()

# Extract carry capacity (k) and growth rate (r) for one species (change 'a' to the desired species)
k <- exp(coef(gam_model)[1 + (which(levels(data$Species) == "a") - 1) * 2])
r <- coef(gam_model)[2 + (which(levels(data$Species) == "a") - 1) * 2]

# Display results for one species (change 'a' to the desired species)
cat("Carry Capacity (k) for Species 'a':", k, "\n")
cat("Growth Rate (r) for Species 'a':", r, "\n")


rm(list = ls())
# Install and load required packages
# install.packages(c("mgcv", "ggplot2"))
library(mgcv)
library(ggplot2)

# Simulated data with 6 species, 5 temperature levels, and 5 nutrients levels
set.seed(123)

# Number of observations per combination
n <- 100

# Create a data frame with 6 species
species <- rep(letters[1:6], each = 5 * n)

# Create a data frame with 5 temperature levels and 5 nutrients levels
temp <- rep(1:5, each = n)
nutrients <- rep(1:5, times = n)

# Simulate response variable for each species
response <- rnorm(length(species), mean = 10 + 2 * temp + 3 * nutrients, sd = 5)

# Create a data frame with numeric variables
data <- data.frame(Species = factor(species),
                   Temp = rep(temp, times = 6),
                   Nutrients = rep(nutrients, times = 6),
                   Response = response)

# Fit GAM with tensor product for the factorial design
gam_model <- gam(Response ~ te(Temp, Nutrients, by = Species), data = data)



Generate Simulated Data
# Simulated data with 6 species, 5 temperature levels, and 5 nutrients levels
set.seed(123)

# Number of observations per combination
n <- 100

# Create a data frame with 6 species
species <- rep(letters[1:6], each = 5 * n)

# Create a data frame with 5 temperature levels and 5 nutrients levels
temp <- rep(1:5, each = n)
nutrients <- rep(1:5, times = n)

# Simulate response variable for each species
response <- rnorm(length(species), mean = 10 + 2 * temp + 3 * nutrients, sd = 5)

# Create a data frame with numeric variables
data <- data.frame(Species = factor(species),
                   Temp = rep(temp, times = 6),
                   Nutrients = rep(nutrients, times = 6),
                   Response = response)

# Fit GAM with tensor product for the factorial design
gam_model <- gam(Response ~ te(Temp, Nutrients, by = Species), data = data)

# Use predict to extract the fitted values
fitted_values <- predict(gam_model, newdata = data)

# Combine fitted values with the original data
results_tibble <- data %>%
  bind_cols(fitted_values = fitted_values) %>%
  mutate(K = exp(fitted_values), r = fitted_values) %>%
  select(Species, K, r)

# Display the tibble
print(results_tibble)








rm(list = ls())
# Install and load required packages
# install.packages(c("mgcv", "ggplot2", "tidyverse"))
library(mgcv)
library(ggplot2)
library(tidyverse)

# Function to simulate logistic growth
simulate_growth <- function(temp, nutrients, days = 10, K_true, r_true) {
  time_points <- seq(1, days, by = 1)
  carrying_capacity <- K_true
  growth_rate <- r_true
  initial_density <- 1
  
  # Simulate logistic growth
  density <- carrying_capacity / (1 + ((carrying_capacity - initial_density) / initial_density) * exp(-growth_rate * time_points))
  
  # Add some random noise
  density <- density + rnorm(length(time_points), sd = 0.1)
  
  # Create a data frame
  data.frame(Temp = rep(temp, each = days),
             Nutrients = rep(nutrients, each = days),
             Day = rep(time_points, times = length(temp)),
             Density = density)
}

# Simulate data for the 5x5 full factorial experiment
set.seed(123)
species <- rep(letters[1:6], each = 5)
temp_levels <- 1:5
nutrient_levels <- 1:5
days <- 10

simulated_data <- expand.grid(Species = species, Temp = temp_levels, Nutrients = nutrient_levels) %>%
  group_by(Species, Temp, Nutrients) %>%
  do(simulate_growth(.$Temp, .$Nutrients, days = days, K_true = rnorm(1, mean = 100, sd = 20), r_true = rnorm(1, mean = 0.1, sd = 0.02)))

# Fit GAM with tensor product for the factorial design
gam_model <- gam(Density ~ te(Temp, Nutrients, by = Species) + s(Day, k = 5), data = simulated_data)

# Use predict to extract the fitted values
fitted_values <- predict(gam_model, newdata = simulated_data, type = "response")

# Combine fitted values with the original data
results_tibble <- simulated_data %>%
  bind_cols(fitted_values = fitted_values) %>%
  group_by(Species) %>%
  summarise(K = max(fitted_values), r = max(diff(log(fitted_values))))

# Display the tibble
print(results_tibble)




# In this code, I've created a simulate_growth function to simulate logistic growth 
# for each combination of temperature, nutrients, and species. The resulting data is 
# then used to fit a GAM model, and carry capacity (K) and growth rate (r) are extracted 
# using tidy functions.

