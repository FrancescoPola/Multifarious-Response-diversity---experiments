# Function to get directional derivatives from GAMS
## x is a community with the following columns: species (chr), pars (list), E1 (dbl), E2 (dbl), and rate (dbl)
get_directionals <- function(x, new_data, refs){ 
  nested_gams <- x %>% 
    nest(cols =-species) %>% 
    mutate(
      gams = map(cols, ~ gam(rate ~  te(E1, E2),
                             data = .x,
                             method = "REML")),
      predicted = map(gams, ~ predict(.x, newdata = new_data))
    )
  
  #list of gams
  m_list <- (nested_gams$gams)
  #list of spp names
  my_spp_names <- (nested_gams$species)
  # get partial derivatives
  (pd_list <- modify_depth(m_list, 1, ~ get_partials(., refs)))
  # from list to tibble
  (pd_spp <- tibble(
    E1_ref = map(pd_list, "E1"),
    E2_ref = map(pd_list, "E2"),
    pd_E1 = map(pd_list, "pd_E1"),
    pd_E2 = map(pd_list, "pd_E2")) %>%
      dplyr::mutate(sp = my_spp_names) %>%
      relocate(sp, E1_ref, E2_ref, pd_E1, pd_E2) %>%
      unnest(E1_ref, E2_ref, pd_E1, pd_E2))
  # add time
  (pd_spp <-cbind(pd_spp, refs$time) %>% 
      dplyr::rename(time = "refs$time"))
  
  # calculation next value for directional derivatives, and get directional derivatives
  (pd_spp <- pd_spp %>% transform(nxt_value_E1 = c(E1_ref[-1], NA)) %>%
      transform(nxt_value_E2 = c(E2_ref[-1], NA)) %>%
      dplyr::mutate(del_E1 = nxt_value_E1 - E1_ref,
                    del_E2 = nxt_value_E2 - E2_ref,
                    unit_vec_mag =  sqrt(del_E1^2 + del_E2^2),
                    uv_E1 = del_E1 / unit_vec_mag,
                    uv_E2 = del_E2 / unit_vec_mag,
                    dir_deriv = pd_E1 * uv_E1 +  pd_E2 * uv_E2) %>% 
      filter(time != max(time)) %>% 
      dplyr::select(sp, time, E1_ref, E2_ref, dir_deriv))
  
}
