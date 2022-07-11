rm(list = ls())

source("main_model_loop.R")

### set bias parameters
p0_offset <- c(-0.1, 0, 0.1)
sigma_mult <- 1/sqrt(c(0.5, 1, 2))

params <- list("p0" = 0.5,
               "numagents" = 100,
               "numouter" = 100,
               "T_pre" = 1,
               "mu_h" = 50,
               "mu_l" = -50,
               "sigma" = 100,
               "alpha" = 0.1,
               "c0" = 153.7,
               "numperiods" = c(400, 400, 20),
               "dT" = c(0.01, 0.1, 1)
               )

perf_out <- matrix(vector(mode = "list", length = 1), nrow = length(p0_offset), ncol = length(sigma_mult))
for (i in 1:length(p0_offset)) {
  params[["p0_offset"]] <- p0_offset[i]
  for (j in 1:length(sigma_mult)) {
    params[["sigma_mult"]] <- sigma_mult[j]
    perf_out[[i, j]] <- main_model_loop(params)
  }
}
