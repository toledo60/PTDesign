Rcpp::sourceCpp(here::here("src", "maximinOofA.cpp"))

position_matrix <- function(x) {
  m <- ncol(x)
  n <- nrow(x)
  y <- matrix(0, nrow = n, ncol = m)

  row_indices <- rep(1:n, times = m)
  col_indices <- as.vector(x) - min(x) + 1
  values <- rep(1:m, each = n)

  y[cbind(row_indices, col_indices)] <- values
  y
}

maximin_oofa <- function(
  n,
  k,
  num_chains = 11,
  total_iterations = 1000000,
  swap_interval = 10,
  alpha = 0.7,
  T0 = 1e-04,
  no_improvement = 5000
) {

  start_time <- Sys.time()
  x <- position_matrix(random_oofa(n,k))

  out <- pt_optim(design = x, M = num_chains, Nmax = total_iterations,
    Nswap = swap_interval, alpha = alpha, T0=T0, tolerance = no_improvement)
  
  end_time <- Sys.time()

  return(list(
    design = out$design,
    measure = out$measure,
    temp_schedule = out$t0,
    elapsed_time = end_time - start_time,
    ntoal = out$ntotal
  ))
}
