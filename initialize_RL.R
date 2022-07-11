initialize_RL <- function(r, mu_h, mu_l, sigma, alpha, p0, numagents, numperiods, dT, p0_offset, sigma_mult, c0) {
  
  if (!("Z" %in% names(r))) {
    if ("fixed_Z" %in% names(r)) {
      r[["Z"]] <- r[["fixed_Z"]]
    } else {
      r[["Z"]] <- matrix(runif(numagents) < p0, nrow = numagents, ncol = 1)
    }
  }
  
  if ("fixed_predraws" %in% names(r)) {
    r[["predraws"]] <- r[["fixed_predraws"]]
  } else {
    r[["predraws"]] <- matrix(rnorm(numagents), nrow = numagents, ncol = 1)
  }
  
  r[["mu_h"]] <- mu_h
  r[["mu_l"]] <- mu_l
  r[["sigma"]] <- sigma
  r[["sigma_est"]] <- sigma*sigma_mult
  r[["alpha"]] <- alpha
  r[["numagents"]] <- numagents
  r[["dT"]] <- dT
  r[["T"]] <- 0
  r[["entrant"]] <- matrix(TRUE, nrow = numagents, ncol = 1)
  r[["P_exit"]] <- matrix(NaN, nrow = numagents, ncol = 1)
  r[["H0_id"]] <- matrix(0, nrow = numagents, ncol = 1)
  r[["L0_id"]] <- matrix(0, nrow = numagents, ncol = 1)
  r[["HE_id"]] <- matrix(0, nrow = numagents, ncol = 1)
  r[["LE_id"]] <- matrix(0, nrow = numagents, ncol = 1)
  r[["HL_id"]] <- matrix(0, nrow = numagents, ncol = 1)
  r[["LL_id"]] <- matrix(0, nrow = numagents, ncol = 1)
  r[["p0"]] <- p0 + p0_offset
  r[["c0"]] <- c0
  r[["P"]] <- matrix(NaN, nrow = numagents, ncol = numperiods + 1)
  r[["P"]][, 1] <- r[["p0"]]
  r[["X"]] <- matrix(NaN, nrow = numagents, ncol = numperiods + 1)
  r[["X"]][, 1] <- 0
  r[["D"]] <- matrix(NaN, nrow = numagents, ncol = numperiods + 1)
  r[["D"]][, 1] <- 0
  r[["G"]] <- matrix(NaN, nrow = numagents, ncol = numperiods + 1)
  r[["G"]][, 1] <- 1
  r["iter"] <- 2
  
  gamma <- sqrt(1 + 8*alpha*(r[["sigma_est"]]^2) / ((mu_h - mu_l)^2))
  
  r[["p"]] <- (-mu_l)*(gamma - 1) / (mu_h*(gamma + 1) + (-mu_l)*(gamma - 1))
  r[["T_exit"]] <- matrix(NaN, nrow = numagents, ncol = 1)
  r[["numHH"]] <- matrix(NaN, nrow = 1, ncol = numperiods + 1)
  r[["numLL"]] <- matrix(NaN, nrow = 1, ncol = numperiods + 1)
  r[["numHL"]] <- matrix(NaN, nrow = 1, ncol = numperiods + 1)
  r[["numLH"]] <- matrix(NaN, nrow = 1, ncol = numperiods + 1)
  r[["numLL"]][, 1] <- 0
  r[["numHL"]][, 1] <- 0
  
  return(r)
  
}