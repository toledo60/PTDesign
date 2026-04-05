maxpro_temps <- function(n, k, M, alpha) {
  crit1 <- 1.0 / (n - 1)
  crit2 <- (1.0 / ((n - 1)^(k - 1) * (n - 2)))^(1.0 / k)
  delta0 <- crit2 - crit1
  T0 <- -delta0 * (1.0 / log(0.99))
  sapply(0:(M - 1), function(i) T0 * alpha^((M - 1) - i))
}

#' Maximum Projection Latin Hypercube Designs
#'
#' Constructs a maximum projection (MaxPro) Latin hypercube design using a
#' parallel tempering optimization algorithm.
#'
#' @param n Integer. Number of design points.
#' @param k Integer. Number of factors.
#' @param num_chains Integer. Number of parallel tempering replicas.
#' @param total_iterations Integer. Maximum number of optimization iterations.
#' @param alpha Numeric. Cooling rate used when constructing a default
#'   temperature schedule.
#' @param swap_interval Integer. Number of iterations between replica swap
#'   attempts.
#' @param no_improvement Integer. Terminate the algorithm after this many
#'   consecutive iterations without improvement.
#' @param temperatures Optional numeric vector of length \code{num_chains}
#'   giving the temperature schedule. If \code{NULL}, a geometric schedule is computed.
#'
#' @return A list with components:
#' \describe{
#'   \item{design}{An \code{n x k} numeric matrix containing the optimized MaxPro
#'   Latin hypercube design.}
#'   \item{t0}{Initial temperature schedule used by the algorithm.}
#'   \item{measure}{The final value of the MaxPro criterion.}
#'   \item{ntotal}{Total number of iterations performed.}
#' }
#'
#' @examples
#' set.seed(123)
#' maxpro_lhd(n = 20, k = 4)
#'
#' @export
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
