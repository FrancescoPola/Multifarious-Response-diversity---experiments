multi_performance_plot <- function(comm){
  these_E2 <- seq( min(comm$E2), max(comm$E2), length.out= 6)
  temp <- comm %>%
    filter(E2 %in% these_E2) %>%
    mutate(E2 = as.factor(E2))
  fig1 <- temp %>% 
    ggplot( aes(x = E1, y = rate, group = as.factor(E2), col = as.factor(E2))) +
    geom_point()+
    geom_smooth(se = FALSE) +
    facet_wrap(~species, ncol = 1, nrow = s) +
    theme_bw()+
    scale_color_viridis_d(end = 0.2,begin = 0.9, option = 'inferno') +
    labs(x = "E1", y = "rate", tag = "(a)")
  
  these_E1 <- seq( min(comm$E1), max(comm$E1), length.out= 6)
  temp <- comm %>%
    filter(E1 %in% these_E1) %>%
    mutate(E1 = as.factor(E1))
  fig2 <- temp %>% 
    ggplot( aes(x = E2, y = rate, group = as.factor(E1), col = as.factor(E1))) +
    geom_point()+
    geom_smooth(se = FALSE) +
    facet_wrap(~species, ncol = 1, nrow = s) +
    theme_bw()+
    scale_color_viridis_d(end = 0.2,begin = 0.9, option = 'inferno') +
    labs(x = "E2", y = "rate", tag = "(b)")
  
  
  performance_comm1 <- fig1 + fig2
  performance_comm1
  
}