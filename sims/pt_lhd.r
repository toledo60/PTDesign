library(dplyr)
library(parallel)
library(PTDesign)

# Setup replication grid

m <- c(5,10,20)
n_multiples <- c(2,5,10)

single_grid <- expand.grid(m=m, n_multiples = n_multiples) |> 
  mutate(n = m*n_multiples) |> 
  arrange(m,n)

alpha <- c(0.25, 0.5, 0.7, 0.9)
n_chains <- c(7, 11, 15)

result <- do.call(
  rbind,
  lapply(1:nrow(single_grid), function(i) {
    cbind(
      single_grid[i, , drop = FALSE],
      expand.grid(alpha = alpha, n_chains = n_chains)
    )
  })
)
reps <- 100
dat <- do.call(rbind, replicate(reps, result, simplify = FALSE))

# maximin LHD Simulation

num_cores <- parallel::detectCores() - 1

measures <- mclapply(
  1:nrow(dat),
  function(i) {
    maximin_lhd(
      n = dat$n[i],
      k = dat$m[i],
      num_chains = dat$n_chains[i],
      p = 15,
      swap_interval = 10,
      total_iterations = 1000000,
      no_improvement = 5000,
      alpha = dat$alpha[i]
    )$measure
  },
  mc.cores = num_cores
)

measures <- unlist(measures)
final_dat <- cbind(dat, measure = measures)

saveRDS(final_dat, file = here::here("sims","data","PT_maximinLHD_sim.rds"))

# MaxPro LHD Simulation

measures <- mclapply(
  1:nrow(dat),
  function(i) {
    maxpro_lhd(
      n = dat$n[i],
      k = dat$m[i],
      num_chains = dat$n_chains[i],
      swap_interval = 10,
      total_iterations = 1000000,
      no_improvement = 5000,
      alpha = dat$alpha[i]
    )$measure
  },
  mc.cores = num_cores
)

measures <- unlist(measures)
final_dat <- cbind(dat, measure = measures)

saveRDS(final_dat, file = here::here("sims","data","PT_maxproLHD_sim.rds"))