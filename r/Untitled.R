
## Modify our function "get_partials" so that you can specify your variable names.


get_partials <- function(m, refs) {
  refs$pd_temperature <- NA
  refs$pd_nutrients <- NA
  for(i in 1:nrow(refs)) {
    refs$pd_temperature[i] <- partial_derivatives(m,
                                                  data = refs[i,],
                                                  type = "central",
                                                  focal = "temperature")$partial_deriv
    refs$pd_nutrients[i] <- partial_derivatives(m,
                                                data = refs[i,],
                                                type = "central",
                                                focal = "nutrients")$partial_deriv
  }
  refs
}