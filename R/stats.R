#' Convert Probability Masses to a Complementary Cumulative Distribution
#'
#' Convert a vector of probability masses into its complementary
#' cumulative distribution function (CCDF).
#'
#' For each element of \code{h}, the corresponding CCDF value is
#' obtained by summing that element and all subsequent elements:
#' \deqn{G_i = \sum_{j=i}^{m} h_j,}
#' where \eqn{m} is the length of \code{h}.
#'
#' @param h A numeric vector of probability masses, typically
#'   representing a discrete probability distribution.
#'
#' @return A numeric vector of the same length as \code{h}.
#'   The \code{i}-th element is the sum of \code{h[i:length(h)]}.
#'
#' @examples
#' h <- c(0.2, 0.3, 0.5)
#' prob2ccdf(h)
#'
#' @export
prob2ccdf <- function(h){
  rev(cumsum(rev(h)))
}


#' Normalize Cumulative Distributions
#'
#' Normalize each column of a matrix by its final cumulative value.
#' This is useful for converting cumulative sums into cumulative
#' distribution functions (CDFs).
#'
#' @param rand A numeric matrix containing cumulative values.
#'
#' @return A numeric matrix with the same dimensions as \code{rand}.
#'   Each column is divided by its value in the last row, so that
#'   the last row of the returned matrix is equal to one.
#'
#' @examples
#' x <- matrix(c(1, 2, 3, 2, 4, 6), nrow = 3)
#' normalize_cdf(x)
#'
#' @export
normalize_cdf <- function(rand){
  sweep(rand, 2, rand[nrow(rand),], "/")
}


#' Generate Random Cumulative Distributions
#'
#' Generate random cumulative distributions from a collection of
#' Gamma-distributed random variables.
#'
#' Random values are generated using \code{\link{rand_gamma}} and
#' accumulated within each generated sample. By default, the
#' resulting values are unnormalized cumulative sums. If
#' \code{normalize} is \code{TRUE}, each cumulative sum is
#' normalized so that its final value is one.
#'
#' @param alpha A parameter passed to \code{\link{rand_gamma}}.
#' @param np An integer specifying the number of random samples
#'   to generate.
#' @param normalize A logical value indicating whether the
#'   cumulative sums should be normalized to form CDFs.
#'   Defaults to \code{FALSE}.
#'
#' @return A numeric matrix containing the cumulative sums of
#'   the generated random values. If \code{normalize} is
#'   \code{TRUE}, each column is normalized to have a final
#'   value of one.
#'
#' @seealso \code{\link{rand_gamma}},
#'   \code{\link{normalize_cdf}}
#'
#' @export
rand_cdf <- function(alpha, np, normalize=FALSE){
  out = rand_gamma(alpha, np=np)
  out = apply(out, 1, cumsum)
  if(normalize){
    out = normalize_cdf(out)
  } 
  return(out)
}


#' Compute the Cumulative Column Marginal
#'
#' Compute the cumulative sum of the column-wise sums of a matrix.
#'
#' @param X A numeric matrix.
#'
#' @return A numeric vector whose \code{i}-th element is the
#'   cumulative sum of the column sums of \code{X} through column
#'   \code{i}.
#'
#' @examples
#' X <- matrix(1:6, nrow = 2)
#' colMarginal(X)
#'
#' @export
colMarginal <- function(X){
  cumsum(colSums(X))  
}


#' Compute Simulation-Based P-values and Confidence Intervals
#'
#' Compute simulation-based p-values and confidence intervals by
#' comparing observed values with values obtained from simulated
#' distributions.
#'
#' For each value in \code{brks}, the corresponding row of simulated
#' values is selected according to the breakpoints in \code{est$breaks}.
#' The simulated values are then compared with the corresponding
#' value in \code{tvalue} to obtain a two-sided p-value.
#'
#' @param rand A numeric matrix containing simulated values. Rows
#'   correspond to breakpoints and columns correspond to simulated
#'   samples.
#' @param est An object containing a numeric vector \code{breaks},
#'   specifying the breakpoints corresponding to the rows of
#'   \code{rand}.
#' @param brks A numeric vector of breakpoints at which p-values
#'   and confidence intervals are to be evaluated.
#' @param tvalue A numeric vector of observed or estimated values
#'   corresponding to the breakpoints in \code{brks}.
#'
#' @return A data frame containing the following columns:
#'   \describe{
#'     \item{x}{The breakpoints supplied in \code{brks}.}
#'     \item{pvalue}{The simulation-based two-sided p-values.}
#'     \item{lower}{The lower quantile of the simulated distribution.}
#'     \item{upper}{The upper quantile of the simulated distribution.}
#'     \item{estimates}{The median of the simulated distribution.}
#'   }
#'
#' @details
#' For each breakpoint, the simulated values are sorted and the
#' empirical position of the corresponding value in \code{tvalue}
#' is used to determine a two-sided p-value. Confidence intervals
#' are obtained from the empirical distribution of the simulated
#' values using probabilities \code{pvalue}, \code{1 - pvalue},
#' and \code{0.5}.
#'
#' @examples
#' \dontrun{
#' result <- simulated_pvalue(
#'   rand = simulated_values,
#'   est = fit,
#'   brks = fit$breaks,
#'   tvalue = fit$estimate
#' )
#' }
#'
#' @export
simulated_pvalue <-function(rand, est, brks, tvalue){
  CI <- matrix(0, nrow = length(brks), ncol = 3)
  colnames(CI) <- c("lower", "upper", "estimates")
  
  #xv = numeric(length(brks))
  pv <- numeric(length = length(brks))
  for(i in 1:length(brks)){
    ind = which.max(est$breaks > brks[i])-1
    ps <- seq(0, 1, length.out=ncol(rand))
    rs <- sort(rand[ind,])
    lower <- ps[which.min(rs < tvalue[i])]
    upper <- rev(ps)[which.min(tvalue[i]<rev(rs))]
    pv[i] <- min(lower, upper)
    
    #xv[i] <- est$breaks[ind]
    ci <- quantile(rand[ind,], probs = c(pv[i], 1-pv[i], 0.5))
    CI[i,] <- ci
  }
  data.frame(x=brks, pvalue=pv, CI)
}
