Rcpp::sourceCpp(here::here("src", "maximin.cpp"))

maximin_lhd <- function(
  n,
  k,
  num_chains = 7,
  total_iterations = 1000000,
  swap_interval = 10,
  alpha = 0.25,
  p = 15,
  no_improvement = 5000
) {

  start_time <- Sys.time()
  out <- maximinLHD(
    n = n,
    k = k,
    M = num_chains,
    Nmax = total_iterations,
    Nswap = swap_interval,
    alpha = alpha,
    p = p,
    tolerance = no_improvement
  )
  end_time <- Sys.time()

  return(list(
    design = out$design,
    measure = out$measure,
    temp_schedule = out$t0,
    elapsed_time = end_time - start_time,
    ntoal = out$ntotal
  ))
}
