
get_optima <- function(x, new_data){ 
  nested_gams <- x %>% 
    nest(cols =-species) %>% 
    mutate(
      gams = map(cols, ~ gam(rate ~  te(E1, E2),
                             data = .x,
                             method = "REML")),
      predicted = map(gams, ~ predict(.x, newdata = new_data))
    )
  
  rates <- nested_gams %>% unnest(cols)
  
  rates %>% group_by(species) %>% 
    filter(rate == max(rate)) %>% 
    dplyr::select(E1, E2)
}
