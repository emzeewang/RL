iterate_RL <- function(r) {
  
  r[["T"]] <- r[["T"]] + r[["dT"]]
  
  num_high <- sum(r[["Z"]])
  num_low <- r[["numagents"]] - num_high
  
  survivors <- r[["G"]][, r[["iter"]]-1] == 1
  
  if ("draws" %in% names(r)) {
    low_draws <- r[["draws"]][1:num_low, r[["iter"]]-1]
    high_draws <- r[["draws"]][(num_low+1):r[["numagents"]], r[["iter"]]-1]
  } else {
    low_draws <- matrix(rnorm(num_low), nrow = num_low, ncol = 1)
    high_draws <- matrix(rnorm(num_high), nrow = num_high, ncol = 1)
  }
  
  r[["D"]][r[["Z"]] == 0, r[["iter"]]] <- r[["mu_l"]]*r[["dT"]] + sqrt(r[["dT"]])*r[["sigma"]]*low_draws
  r[["D"]][r[["Z"]] == 1, r[["iter"]]] <- r[["mu_h"]]*r[["dT"]] + sqrt(r[["dT"]])*r[["sigma"]]*high_draws
  r[["D"]][!survivors, r[["iter"]]] <- 0
  r[["X"]][, r[["iter"]]] <- r[["X"]][, r[["iter"]]-1] + r[["D"]][, r[["iter"]]]
  r[["P"]][, r[["iter"]]] <- (1 + ((1-r[["p0"]])/r[["p0"]])*exp((-(r[["mu_h"]]-r[["mu_l"]]) / (r[["sigma_est"]]^2)) * (r[["X"]][, r[["iter"]]]-(r[["mu_h"]]+r[["mu_l"]])*T/2)))^(-1)
  
  r[["G"]][r[["G"]][, r[["iter"]]-1] == 0, r[["iter"]]] <- 0
  r[["G"]][r[["G"]][, r[["iter"]]-1] == 1, r[["iter"]]] <- (r[["P"]][r[["G"]][, r[["iter"]]-1] == 1, r[["iter"]]] > r[["p"]])
  
  exit <- (r[["G"]][, r[["iter"]]] == 0) & (r[["G"]][, r[["iter"]]-1] == 1)
  
  r[["T_exit"]][exit, 1] <- r[["T"]]
  r[["P_exit"]][exit, 1] <- r[["P"]][exit, r[["iter"]]]
  
  r[["numHH"]][1, r[["iter"]]] <- sum((r[["Z"]] == 1) & (r[["G"]][, r[["iter"]]] == 1) & (r[["entrant"]]))
  r[["numLH"]][1, r[["iter"]]] <- sum((r[["Z"]] == 0) & (r[["G"]][, r[["iter"]]] == 1) & (r[["entrant"]]))
  r[["numLL"]][1, r[["iter"]]] <- sum((r[["Z"]] == 0) & (r[["G"]][, r[["iter"]]] == 0) & (r[["entrant"]]))
  r[["numHL"]][1, r[["iter"]]] <- sum((r[["Z"]] == 1) & (r[["G"]][, r[["iter"]]] == 0) & (r[["entrant"]]))
  
  r[["iter"]] <- r[["iter"]] + 1
  
  return(r)
  
}