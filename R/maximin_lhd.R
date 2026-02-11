#' Maximin Latin Hypercube Designs
#'
#' Constructs a maximin Latin hypercube design using a parallel tempering
#' optimization algorithm.
#'
#' @param n Integer. Number of design points.
#' @param k Integer. Number of factors.
#' @param num_chains Integer. Number of parallel tempering replicas.
#' @param total_iterations Integer. Maximum number of iterations.
#' @param swap_interval Integer. Number of iterations between replica swap attempts.
#' @param alpha Numeric. Cooling rate controlling the temperature schedule.
#' @param p Integer. Power parameter in the average reciprocal inter-point
#'   distance criterion. Default is 15.
#' @param no_improvement Integer. Terminate the algorithm after this many
#'   consecutive iterations without improvement.
#'
#' @return A list with components:
#' \describe{
#'   \item{design}{An \code{n x k} numeric matrix containing the optimized
#'   Latin hypercube design.}
#'   \item{measure}{The final value of the maximin criterion.}
#'   \item{temp_schedule}{Initial temperature schedule used by the algorithm.}
#'   \item{elapsed_time}{Total elapsed computation time.}
#'   \item{ntoal}{Total number of iterations performed.}
#' }
#'
#' @examples
#' set.seed(123)
#' maximin_lhd(n = 20, k = 3)
#'
#' @export
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
