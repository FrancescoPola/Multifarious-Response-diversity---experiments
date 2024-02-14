#' Plot performance curves, made by `Make_expt()`
#'
#' @param expt A data frame with three columns: E1, E2, and rate. Returned by the `Make_expt()` function.
#' @return patchwork object
#' @export
Plot_performance_curves <- function(expt) {
  
  these_E2 <- expt$E2[seq(1, length(expt$E2), length=5)]
  temp <- expt %>%
    filter(E2 %in% these_E2) %>%
    mutate(E2 = as.factor(E2))
  fig1 <- ggplot(temp, aes(x = E1, y = rate, col = E2)) +
    geom_line()   
  
  these_E1 <- expt$E1[seq(1, length(expt$E2), length=5)]
  temp <- expt %>%
    filter(E1 %in% these_E1) %>%
    mutate(E1 = as.factor(E1))
  fig2 <- ggplot(temp, aes(x = E2, y = rate, col = E1)) +
    geom_line()
  
  list(fig1, fig2)
  
}
