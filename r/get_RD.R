
### Calculating  Response Diversity for all community compositions 
# requires specifying your variables names, example:
# E1_var <- "temperature"
# E2_var <- "nutrients"
# Create a function to perform the calculations
get_RD <- function(species_names, nested_gams, E1_var, E2_var) {
  species_data <- nested_gams %>% filter(species %in% species_names)
  
  m_list <- species_data$gams
  my_spp_names <- paste0("s_", species_data$species)
 
  # Extract partial derivatives
  pd_list <- modify_depth(m_list, 1, ~ get_partials_mod(., env_data, E1_var, E2_var))
  # from list to tibble
  (pd_spp <- tibble(
    E1_ref = map(pd_list, E1_var),
    E2_ref = map(pd_list, E2_var),
    pd_E1 = map(pd_list, paste0("pd_", E1_var)),
    pd_E2 = map(pd_list, paste0("pd_", E2_var)) 
  )%>% 
      dplyr::mutate(sp = my_spp_names) %>%
      unnest(c(E1_ref, E2_ref, pd_E1, pd_E2)))
  
  # add time
  (pd_spp <-cbind(pd_spp, env_data$time) %>% 
      dplyr::rename(time = "env_data$time"))
  
  # calculation next value for directional derivatives, and get directional derivatives
  (pd_spp <- pd_spp %>% transform( nxt_value_E1 = c(E1_ref[-1], NA)) %>%
      transform(nxt_value_E2 = c(E2_ref[-1], NA)) %>%
      dplyr::mutate(del_E1 = nxt_value_E1 - E1_ref,
                    del_E2 = nxt_value_E2 - E2_ref,
                    unit_vec_mag =  sqrt(del_E1^2 + del_E2^2),
                    uv_E1 = del_E1 / unit_vec_mag,
                    uv_E2 = del_E2 / unit_vec_mag,
                    dir_deriv = pd_E1 * uv_E1 +  pd_E2 * uv_E2) %>% 
      filter(time != max(time)))
  
  
  red_spp <- pd_spp %>%  dplyr::select(sp, time, E1_ref, E2_ref, dir_deriv)
  # from long to wide
  rdiv_1 <- red_spp %>%
    spread( sp, dir_deriv)
  
  rdiv_1[is.na(rdiv_1)] <- 0
  
  
  # actual calculation for only the same species used above
  rdiv_1$dissimilarity<-apply(dplyr:: select(rdiv_1, starts_with("s")), 1, resp_div, sign_sens = F)
  rdiv_1$divergence<-apply(dplyr:: select(rdiv_1, starts_with("s")), 1, resp_div, sign_sens = T)                  
  # Add species names to the result
  rdiv_1$species <- paste(species_names, collapse = "_")
  rdiv_1 <- rdiv_1 %>% select(time, E1_ref, E2_ref, dissimilarity, divergence, species)
  # community 1 
  # Return the result
  return(rdiv_1)
}

