compute_continuous_discount <- function(dT, numperiods, alpha) {
  
  tvec <- seq(0, dT[1]*(numperiods[1]-1), dT[1])
  avc <- exp(-alpha*tvec) * ((1 - exp(-alpha*dT[1])) / alpha / dT[1])
  
  avec <- c(1, avc)
  
  for (i in 2:length(dT)) {
    tvec <- seq(tail(tvec, 1) + dT[i-1], tail(tvec, 1) + dT[i-1] + dT[i]*(numperiods[i]-1), dT[i])
    avc <- exp(-alpha*tvec) * ((1 - exp(-alpha*dT[i])) / alpha / dT[i])
  
    avec <- c(avec, avc)
  }
  
  return(avec)
  
}