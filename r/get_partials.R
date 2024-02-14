
## Modify our function "get_partials" so that you can specify your variable names.
# change E1_var = "E1", E2_var = "E2", with  E1_var = "your variable 1 name" & E2_var = "your variable 2 name"

get_partials <- function(m, refs, E1_var = "E1", E2_var = "E2") {
  for (i in 1:nrow(refs)) {
    refs[[paste0("pd_", E1_var)]][i] <- partial_derivatives(m,
                                                            data = refs[i,],
                                                            type = "central",
                                                            focal = E1_var)$partial_deriv
    refs[[paste0("pd_", E2_var)]][i] <- partial_derivatives(m,
                                                            data = refs[i,],
                                                            type = "central",
                                                            focal = E2_var)$partial_deriv
  }
  refs
}


