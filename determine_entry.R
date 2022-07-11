determine_entry <- function(r, dT) {
  
  T <- dT
  
  num_high <- sum(r[["Z"]])
  num_low <- r[["numagents"]] - num_high
  
  X <- matrix(NaN, nrow = r[["numagents"]], ncol = 1)
  
  if ("predraws" %in% names(r)) {
    low_draws <- r[["predraws"]][1:num_low, 1]
    high_draws <- r[["predraws"]][(num_low+1):r[["numagents"]], 1]
  } else {
    low_draws <- matrix(rnorm(numagents), nrow = numagents, ncol = 1)
    high_draws <- matrix(rnorm(numagents), nrow = numagents, ncol = 1)
  }
  
  sigma <- r[["sigma"]]
  sigma_est <- r[["sigma_est"]]
  
  X[r[["Z"]] == 0] <- r[["mu_l"]]*dT + sqrt(dT)*sigma*low_draws
  X[r[["Z"]] == 1] <- r[["mu_h"]]*dT + sqrt(dT)*sigma*high_draws
  
  r[["p0"]] <- (1 + ((1-r[["p0"]])/r[["p0"]])*exp((-(r[["mu_h"]]-r[["mu_l"]]) / (sigma_est^2)) * (X-(r[["mu_h"]]+r[["mu_l"]])*T/2)))^(-1)
  r[["P"]][, 1] <- r[["p0"]]

  gamma <- sqrt(1 + 8*r[["alpha"]]*(r[["sigma_est"]]^2) / ((r[["mu_h"]] - r[["mu_l"]])^2))
  
  orig_rp <- r[["p"]]
  r[["p"]] <- (-r[["mu_l"]])*(gamma - 1) / (r[["mu_h"]]*(gamma + 1) + (-r[["mu_l"]])*(gamma - 1))
  
  E_pmu <- r[["p"]]*r[["mu_h"]] + (1-r[["p"]])*r[["mu_l"]]
  E_p0mu <- r[["p0"]]*r[["mu_h"]] + (1-r[["p0"]])*r[["mu_l"]]
  Rp <- (E_p0mu/r[["alpha"]]) - (E_pmu/r[["alpha"]])*((1-r[["p0"]])/(1-r[["p"]]))*((1-r[["p0"]])/(1-r[["p"]])*r[["p"]]/r[["p0"]])^((gamma-1)/2)
  
  r[["entrant"]] <- (Rp > r[["c0"]]) & (r[["p0"]] >= r[["p"]])
  r[["G"]][, 1] <- r[["entrant"]]
  
  r[["numHH"]][, 1] <- sum((r[["Z"]] == 1) & r[["entrant"]])
  r[["numLH"]][, 1] <- sum((r[["Z"]] == 0) & r[["entrant"]])
  
  r[["H_id"]] <- r[["Z"]]
  r[["L_id"]] <- !(r[["Z"]])
  r[["H0_id"]] <- matrix(r[["H_id"]] * !(r[["entrant"]]), nrow = r[["numagents"]], ncol = 1)
  r[["L0_id"]] <- matrix(r[["L_id"]] * !(r[["entrant"]]), nrow = r[["numagents"]], ncol = 1)
  r[["HE_id"]] <- matrix(r[["H_id"]] * r[["entrant"]], nrow = r[["numagents"]], ncol = 1)
  r[["LE_id"]] <- matrix(r[["L_id"]] * r[["entrant"]], nrow = r[["numagents"]], ncol = 1) 
  r[["HL_id"]] <- matrix(0, nrow = r[["numagents"]], ncol = 1) 
  r[["LL_id"]] <- matrix(0, nrow = r[["numagents"]], ncol = 1) 
  
  r[["p"]] <- orig_rp
  
  return(r)
  
}