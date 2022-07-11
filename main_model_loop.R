main_model_loop <- function(params) {
  
  source("compute_continuous_discount.R")
  source("initialize_RL.R")
  source("determine_entry.R")
  source("iterate_RL.R")
  
  ### parameters
  numagents <- params[["numagents"]]
  numouter <- params[["numouter"]]
  mu_l <- params[["mu_l"]]
  mu_h <- params[["mu_h"]]
  p0 <- params[["p0"]]
  p0_offset <- params[["p0_offset"]]
  sigma <- params[["sigma"]]
  T_pre <- params[["T_pre"]]
  alpha <- params[["alpha"]]
  c0 <- params[["c0"]]
  sigma_mult <- params[["sigma_mult"]]
  numperiods <- params[["numperiods"]]
  dT <- params[["dT"]]
  
  ### calculate post-entry time increment variables
  dT1 <- dT[1]
  avec <- compute_continuous_discount(dT, numperiods, alpha)
  numperiods_tot <- numperiods[1]
  
  for (i in 2:length(dT)) {
    numperiods_tot <- numperiods_tot + numperiods[i]
  }
  
  tvec <- 0
  dtvec <- 0
  t <- 2
  
  for (k in 1:length(dT)) {
    for (i in 1:numperiods[k]) {
      tvec[t] <- tvec[t-1] + dT[k]
      dtvec[t] <- dT[k]
      t <- t + 1
    }
  }
  
  ### initialize variables to store outputs
  nHH <- matrix(0, nrow = 1, ncol = numperiods_tot + 1)
  nLL <- matrix(0, nrow = 1, ncol = numperiods_tot + 1)
  nHL <- matrix(0, nrow = 1, ncol = numperiods_tot + 1)
  nLH <- matrix(0, nrow = 1, ncol = numperiods_tot + 1)
  nH0 <- 0
  nL0 <- 0
  
  n_exit <- matrix(0, nrow = 1, ncol = numperiods_tot)
  n_high_exit <- matrix(0, nrow = 1, ncol = numperiods_tot)
  n_low_exit <- matrix(0, nrow = 1, ncol = numperiods_tot)
  
  profit_cum <- 0
  profit_store <- matrix(0, nrow = numagents, ncol = numouter)
  entrant_store <- matrix(0, nrow = numagents, ncol = numouter)
  
  ### main loop
  
  for (m in 1:numouter) {
    
    print(m)
    
    r <- list("init" = NaN)
    r <- initialize_RL(r, mu_h, mu_l, sigma, alpha, p0, numagents, numperiods_tot, dT1, p0_offset, sigma_mult, c0)
    r <- determine_entry(r, T_pre)
    
    for (k in 1:length(dT)) {
      r[["dT"]] <- dT[k]
      for (p in 1:numperiods[k]) {
        r <- iterate_RL(r)
      }
    }
    
    nHH <- nHH + r[["numHH"]]
    nLL <- nLL + r[["numLL"]]
    nHL <- nHL + r[["numHL"]]
    nLH <- nLH + r[["numLH"]]
    nH0 <- nH0 + sum((!r[["entrant"]]) & (r[["Z"]]))
    nL0 <- nL0 + sum((!r[["entrant"]]) & (!r[["Z"]]))
    
    exit <- r[["G"]][, 1:(dim(r[["G"]])[2]-1)] - r[["G"]][, 2:(dim(r[["G"]])[2])]
    n_exit <- n_exit + colSums(exit)
    n_high_exit <- n_high_exit + colSums(exit[r[["Z"]],])
    n_low_exit <- n_low_exit + colSums(exit[!r[["Z"]],])
    
    profit_cum <- profit_cum + sum(r[["D"]] %*% avec)
    profit_store[, m] <- r[["D"]] %*% avec
    entrant_store[, m] <- r[["entrant"]]
    
  }
  
  s <- vector(mode = "list", length = 0)
  s[["nHH"]] <- nHH
  s[["nLL"]] <- nLL
  s[["nHL"]] <- nHL
  s[["nLH"]] <- nLH
  s[["nH0"]] <- nH0
  s[["nL0"]] <- nL0
  s[["profit_cum"]] <- profit_cum
  s[["profit_store"]] <- profit_store
  s[["entrant_store"]] <- entrant_store
  s[["n_exit"]] <- n_exit
  s[["n_high_exit"]] <- n_high_exit
  s[["n_low_exit"]] <- n_low_exit
  s[["profit_avg_entrants"]] <- profit_cum / (nHH[length(nHH)] + nLH[length(nLH)] + nLL[length(nLL)] + nHL[length(nHL)]) - c0
  s[["profit_avg_all"]] <- s[["profit_avg_entrants"]] * (numouter*numagents - nL0 - nH0) / (numouter*numagents)

  return(s)
  
}