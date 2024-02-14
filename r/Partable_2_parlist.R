#' Convert a table of parameters to a list of parameters
#'
#' @param E1_series The series of values of $E_1$ used in the experiment.
#' @param E2_series The series of values of $E_2$ used in the experiment.
#' @param pars The parameter list for the Eppley function
#' @return A data frame with three columns: E1, E2, and rate.
#' @export
Partable_2_parlist <- function(par_table) {
  par_list <- list()
  for(i in 1:nrow(par_table)) 
    par_list[[i]] <- as.list(t(par_table[i,])[,1])
  par_list
}