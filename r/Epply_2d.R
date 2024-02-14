#' Epply functional response curve with 2 environmental factors, with potential interaction
#' effect between the two environmental factors.
#' 
#' Parameters are given in a list, and are:
#' $z$ controls location of maximum.
#' $w$ controls range of $E$ over which the rate is positive.
#' $a$ scaling constant.
#' $b$ controls rate of increase towards the maximum rate, as $E$ increases.
#' $zint$ interaction effect on location of maximum
#' $sd_rate$ is the standard deviation of the noise that is added to the rate
#'
#' @param E1 Value of environmental factor 1
#' @param E2 Value of environmental factor 2
#' @param pars A list of parameters
#' @return The calculated rate
#' @export
Eppley_2d <- function(E1, E2, pars) {
  with(pars,
       {
         rate <- a1 * exp(b1 * E1) *
                         (1 - (((E1 - z1) + zint*(E1 - z1)*(E2 - z2)) / (w1/2))^2) +
                 a2 * exp(b2 * E2) *
                         (1 - ((E2 - z2) / (w2/2))^2) 
         rate
       }
  )
}
