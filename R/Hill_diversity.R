#' Hill Number Diversity
#'
#' Compute the Hill number of order \code{q} for a vector or matrix of counts.
#'
#' For a single community (numeric vector \code{counts}), 
#' \deqn{
#'   ^qD = \begin{cases}
#'     \bigl(\sum_i p_i^q\bigr)^{1/(1-q)}, & q \neq 1,\\
#'     \exp\!\bigl(-\sum_i p_i \ln p_i\bigr),   & q = 1,
#'   \end{cases}
#' }
#' where \eqn{p_i = \text{counts}_i / \sum_i \text{counts}_i}.
#'
#' If \code{counts} is a numeric matrix or data.frame, each row is treated
#' as an independent community and the function returns a numeric vector
#' of the same length as \code{nrow(counts)}.
#'
#' @param counts Numeric vector of non‐negative integer counts, or a
#'   numeric matrix/data.frame where rows are communities and columns species.
#' @param q Non‐negative numeric order of the Hill number.  \code{q=0}
#'   yields species richness, \code{q=1} yields the exponential of Shannon,
#'   \code{q=2} yields the inverse Simpson, etc.
#' @return A numeric value (if \code{counts} is a vector) or numeric vector
#'   (if \code{counts} is a matrix) of Hill diversities.
#' @examples
#' # single vector
#' counts <- c(10, 5, 1, 0, 0)
#' Hill_diversity(counts, q = 0)  # 3 species
#' Hill_diversity(counts, q = 1)  # exp(Shannon)
#' Hill_diversity(counts, q = 2)  # inverse Simpson
#'
#' # multiple communities
#' comms <- matrix(c(10,5,1,  3,3,4), ncol=2, byrow=FALSE)
#' Hill_diversity(comms, q = 1)
#' @export
Hill_diversity <- function(counts, q) {
  if (any(counts < 0)) stop("counts must be non‐negative")
  # handle vector vs matrix
  if (is.vector(counts)) {
    x <- counts
    tot <- sum(x)
    if (tot == 0) return(0)
    p <- x / tot
    if (q == 1) {
      p <- p[p > 0]
      return(exp(-sum(p * log(p))))
    } else if (q == 0) {
      return(sum(x > 0))
    } else {
      return((sum(p^q))^(1 / (1 - q)))
    }
  } else {
    # assume counts is a matrix of count data: columns = communities
    mat <- as.matrix(counts)
    apply(mat, 2, function(x) Hill_diversity(x, q))
  }
}
