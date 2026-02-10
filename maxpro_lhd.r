Rcpp::sourceCpp(here::here("src", "maxpro.cpp"))

maxpro_temps <- function(n, k, M, alpha) {
  crit1 <- 1.0 / (n - 1)
  crit2 <- (1.0 / ((n - 1)^(k - 1) * (n - 2)))^(1.0 / k)
  delta0 <- crit2 - crit1
  T0 <- -delta0 * (1.0 / log(0.99))
  return(sapply(0:(M - 1), function(i) T0 * alpha^(M - i)))
}

maxpro_lhd <- function(
  n,
  k,
  num_chains = 7,
  total_iterations = 1000000,
  alpha = 0.9,
  swap_interval = 10,
  no_improvement = 5000,
  temperatures = NULL
) {
  if (is.null(temperatures)) {
    temperatures <- maxpro_temps(n = n, k = k, M = num_chains, alpha = alpha)
  }
  out <- MaxProLHD(
    n = n,
    k = k,
    M = num_chains,
    Nswap = swap_interval,
    Nmax = total_iterations,
    tolerance = no_improvement,
    temp_sched = temperatures
  )
  out
}
