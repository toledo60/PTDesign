library(dplyr)
library(parallel)

library(SLHD)
library(MaxPro)

# Setup replication grid

m <- c(5,10,20)
n_multiples <- c(2,5,10)

single_grid <- expand.grid(m=m, n_multiples = n_multiples) |> 
  mutate(n = m*n_multiples) |> 
  arrange(m,n)

reps <- 100
dat <- do.call(rbind, replicate(reps, single_grid , simplify = FALSE))

# SLHD: maximin LHD Simulation

num_cores <- parallel::detectCores() - 1

measures <- mclapply(
  1:nrow(dat),
  function(i) {
    SLHD::maximinSLHD(t=1, m= dat$n[i], k = dat$m[i])$measure
  },
  mc.cores = num_cores
)

measures <- unlist(measures)
final_dat <- cbind(dat, measure = measures)

saveRDS(final_dat, file = here::here("sims","data","SLHD_maximin_sim.rds"))

# MaxPro LHD Simulation

num_cores <- parallel::detectCores() - 1

measures <- mclapply(
  1:nrow(dat),
  function(i) {
    MaxPro::MaxProLHD(n = dat$n[i], p = dat$m[i])$measure
  },
  mc.cores = num_cores
)

measures <- unlist(measures)
final_dat <- cbind(dat, measure = measures)

saveRDS(final_dat, file = here::here("sims","data","maxproLHD_sim.rds"))
