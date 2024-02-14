### Calculating Potential Response Diversity for all community compositions 
# requires specifying your variables names, example:
# E1_var <- "temperature"
# E2_var <- "nutrients"
# Create a function to perform the calculations
get_potential_RD <- function(species_names, nested_gams, E1_var, E2_var) {
  
  species_data <- nested_gams %>% filter(species %in% species_names)
  
  m_list <- species_data$gams
  my_spp_names <- species_data$species
 
  # Extract partial derivatives
  pd_list <- modify_depth(m_list, 1, ~ get_partials_mod(., refs2, E1_var, E2_var))
  

  # Convert to tibble
  pd_spp <- tibble(
    E1_ref = map(pd_list, E1_var),
    E2_ref = map(pd_list, E2_var),
    pd_E1 = map(pd_list, paste0("pd_", E1_var)),
    pd_E2 = map(pd_list, paste0("pd_", E2_var))
  ) %>%
    dplyr::mutate(sp = my_spp_names) %>%
    unnest(c(E1_ref, E2_ref, pd_E1, pd_E2))
  
  
  # Calculate diversity
  radius <- 1
  num_arrows <- 100
  pd_spp <- crossing(angle = rep(seq(0, 2*pi, length = num_arrows)),
                     E1_ref = pd_spp$E1_ref,
                     E2_ref = pd_spp$E2_ref) %>%
    full_join(pd_spp) %>%
    mutate(E1 = cos(angle) * radius,
           E2 = sin(angle) * radius,
           dir_deriv = E1 * pd_E1 + E2 * pd_E2,
           unit_vec_mag = sqrt(E1^2 + E2^2))
  
  # Filter unnecessary columns
  pd_spp <- pd_spp %>%
    dplyr::select(angle, sp, E1_ref, E2_ref, dir_deriv)
  
  # Calculate diversity measures
  Div_loc_dir <- pd_spp %>% 
    dplyr::group_by(E1_ref, E2_ref, angle) %>% 
    summarise(div = resp_div(dir_deriv, sign = FALSE))
  
  RDiv <- tibble(Div_loc_dir %>% dplyr::group_by(E1_ref, E2_ref) %>% 
                   summarise(mean = mean(div)) %>% 
                   ungroup() %>% 
                   summarise(dissimilarity = mean(mean)))
  
  Div_loc_dir2 <- pd_spp %>% 
    dplyr::group_by(E1_ref, E2_ref, angle) %>% 
    summarise(div = resp_div(dir_deriv, sign = TRUE))
  
  RDiv1 <- tibble(Div_loc_dir2 %>% dplyr::group_by(E1_ref, E2_ref) %>% 
                    summarise(mean = mean(div)) %>% 
                    ungroup() %>% 
                    summarise(divergence = mean(mean))) %>% 
    cbind(RDiv)
  
  # Add species names to the result
  RDiv1$species <- paste(species_names, collapse = "_")
  
  # Return the result
  return(RDiv1)
}
