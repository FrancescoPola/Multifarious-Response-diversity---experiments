#' Calculate rates for a combination of values of $E_1$ and $E_2$, for a given parameter set
#' Somewhat similar to running an experiment to measure the rates at
#' different combinations of the two environmental factors.
#'
#' @param E1_series The series of values of $E_1$ used in the experiment.
#' @param E2_series The series of values of $E_2$ used in the experiment.
#' @param pars The parameter list for the Eppley function
#' @return A data frame with three columns: E1, E2, and rate.
#' @export
Make_expt <- function(E1_series, E2_series, pars = pars, perf_func = Eppley_2d_v2)
{
  expt <- crossing(E1 = E1_series,
                   E2 = E2_series) %>%
    dplyr::mutate(temp_rate = perf_func(E1, E2, pars),
           noise = rnorm(length(temp_rate), 0, pars$sd_rate),
           rate = temp_rate + noise) %>%
    dplyr::select(-temp_rate, -noise)
  list(expt)
}