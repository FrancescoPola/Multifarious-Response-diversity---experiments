

directional_deriv <- function(E1_p, E2_p, partial_deriv_E1, partial_deriv_E2){
directional_derivative =  (partial_deriv_E1 * (E1_p/(sqrt(E1_p^2 + E2_p^2)))) + ((partial_deriv_E2) * (E2_p/(sqrt(E1_p^2 + E2_p^2))))
}
