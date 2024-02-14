## Function to get community 
get_community <- function(s, species_pars){
  
  par_table <- with(species_pars, {
      tibble(a1 = rnorm(s, a1_mean, 0),#a1_mean/100),
                           b1 = rnorm(s, b1_mean, 0),#b1_mean/100),
                           z1 = runif(s, min = z1_mean - z1_range/2, max = z1_mean + z1_range / 2),
                           w1 = rnorm(s, w1_mean, 0),#w1_mean/100),
                           a2 = rnorm(s, a2_mean, 0),#a2_mean/100),
                           b2 = rnorm(s, b2_mean, 0),#b2_mean/100),
                           z2 = runif(s, min = z2_mean - z2_range/2, max = z2_mean + z2_range / 2),
                           w2 = rnorm(s, w2_mean, 0),#w2_mean/100),
                           z_int21 = rnorm(s, zint_mean, zint_mean/10),
                           sd_rate = rnorm(s, sd_rate_mean, sd_rate_mean/10)
       )
  }
  )
   
  par_list <- Partable_2_parlist(par_table)
  ## add performance curves
  species_pars <- tibble(species = paste0("s", 1:s), pars = par_list) %>%
    rowwise() %>%
    mutate(expt = Make_expt(E1_series, E2_series, pars))
  
  (species_pars <- species_pars %>% unnest(expt))
}