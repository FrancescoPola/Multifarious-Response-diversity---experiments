
# Function to get partial derivatives from GAMS

get_partials <- function(m, refs) {
  refs$pd_E1 <- NA
  refs$pd_E2 <- NA
  for(i in 1:nrow(refs)) {
    refs$pd_E1[i] <- partial_derivatives(m,
                                         data = refs[i,],
                                         type = "central",
                                         focal = "E1")$partial_deriv
    refs$pd_E2[i] <- partial_derivatives(m,
                                         data = refs[i,],
                                         type = "central",
                                         focal = "E2")$partial_deriv
  }
  refs
}