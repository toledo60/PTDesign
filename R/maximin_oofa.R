#' Position Matrix for Order-of-Addition Designs
#'
#' Constructs a position matrix representation of an order-of-addition design.
#' Each row corresponds to a run, and each column indicates the position at
#' which a component appears in the sequence.
#'
#' @param x An integer matrix of dimension \code{n x m} representing an
#'   order-of-addition design, where each row is a permutation of the
#'   components.
#'
#' @return An integer matrix of dimension \code{n x m}. Entry \code{(i, j)}
#'   gives the position at which component \code{j} appears in run \code{i}.
#'
#' @examples
#' x <- matrix(c(1, 2, 3,
#'               2, 3, 1,
#'               3, 1, 2),
#'             nrow = 3, byrow = TRUE)
#' position_matrix(x)
#'
#' @export
position_matrix <- function(x) {
  m <- ncol(x)
  n <- nrow(x)
  y <- matrix(0, nrow = n, ncol = m)

  row_indices <- rep(seq_len(n), times = m)
  col_indices <- as.vector(x) - min(x) + 1
  values <- rep(1:m, each = n)

  y[cbind(row_indices, col_indices)] <- values
  y
}


#' Maximin Order-of-Addition Designs
#'
#' Constructs a maximin order-of-addition (OofA) design using a parallel
#' tempering optimization algorithm.
#'
#' @param n Integer. Number of experimental runs.
#' @param k Integer. Number of components to be ordered.
#' @param num_chains Integer. Number of parallel tempering replicas.
#' @param total_iterations Integer. Maximum number of optimization iterations.
#' @param swap_interval Integer. Number of iterations between replica swap
#'   attempts.
#' @param alpha Numeric. Cooling rate controlling the temperature schedule.
#' @param T0 Numeric. Initial temperature for the parallel tempering algorithm.
#' @param no_improvement Integer. Terminate the algorithm after this many
#'   consecutive iterations without improvement.
#'
#' @return A list with components:
#' \describe{
#'   \item{design}{An \code{n x k} integer matrix representing the optimized
#'   order-of-addition design in position-matrix form.}
#'   \item{measure}{The final value of the maximin objective function.}
#'   \item{temp_schedule}{Initial temperature schedule used by the algorithm.}
#'   \item{elapsed_time}{Total elapsed computation time}
#'   \item{ntoal}{Total number of iterations performed.}
#' }
#'
#' @examples
#' set.seed(123)
#' d <- maximin_oofa(n = 12, k = 4)
#'
#' @export
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
