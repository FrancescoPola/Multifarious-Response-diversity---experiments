


partial_deriv_E1_E2 <- function(E1_p, E2_p, my_gams_list, my_sp_names)
{
  # data slice through te(x,z) holding z == 0.4
  ds1 <- data.frame(E1 = c(E1_p), E2 = E2_p)
  
  # evaluate te(E1,E2) at values of E1 & E2
  sm_E1 <- modify_depth(my_gams_list, 1, ~ smooth_estimates(., smooth = "te(E1,E2)", data = ds1) |>
                          add_confint())
  
  # partial derivatives E1 
  (pd_E1_list <- modify_depth(my_gams_list, 1, ~ partial_derivatives(., data = ds1, type = "central", focal = "E1")))
  
  # from list to tibble
  (pd_E1 <- tibble(
    E1= map(pd_E1_list, "data"),
    partial_deriv = map(pd_E1_list, "partial_deriv")) %>% 
      dplyr::mutate(sp = my_sp_names) %>% 
      relocate(sp, E1, partial_deriv) %>% 
      unnest(sp,E1, partial_deriv) %>% 
      dplyr::rename(partial_deriv_E1 = "partial_deriv") )
  
  
  # data slice through te(x,z) holding z == 0.4
  ds2 <- data.frame(E2 = c(E2_p), E1 = E1_p)
  # evaluate te(E1,E2) at values of E1 & E2
  sm_E2 <- modify_depth(my_gams_list, 1, ~ smooth_estimates(., smooth = "te(E1,E2)", data = ds2) |>
                          add_confint())
  
  # partial derivatives E1 
  (pd_E2_list <- modify_depth(my_gams_list, 1, ~ partial_derivatives(., data = ds2, type = "central", focal = "E2")))
  
  # from list to tibble
  (pd_E2 <- tibble(
    E2= map(pd_E2_list, "data"),
    partial_deriv = map(pd_E2_list, "partial_deriv")) %>% 
      dplyr::mutate(sp = my_sp_names) %>% 
      relocate(sp, E2, partial_deriv) %>% 
      unnest(sp, E2, partial_deriv) %>% 
      dplyr::rename(partial_deriv_E2 = "partial_deriv"))
  
  (pd <- merge(pd_E1, pd_E2, by = "sp") %>% 
      relocate(sp, E1, E2, partial_deriv_E1, partial_deriv_E2))
  
}
