#' Remove the Last Element
#'
#' Remove the last element of a vector.
#'
#' @param x A vector.
#'
#' @return A vector with its last element removed.
#'
#' @keywords internal
rmlast <- function(x)x[-length(x)]


#' Construct a Sorted Vector of Breakpoints
#'
#' Construct a sorted vector of unique breakpoints from the supplied
#' left and right endpoints. An optional value \code{tau} can also
#' be included in the set of breakpoints.
#'
#' @param LE A vector of left endpoints.
#' @param RE A vector of right endpoints.
#' @param LS A vector of left endpoints for the second set.
#' @param RS A vector of right endpoints for the second set.
#' @param tau An optional vector of additional breakpoints.
#'   Defaults to \code{NULL}.
#'
#' @return A sorted vector of unique breakpoints obtained by combining
#'   \code{LE}, \code{RE}, \code{LS}, and \code{RS}, and optionally
#'   \code{tau}.
#'
#' @keywords internal
setbreaks <- function(LE, RE, LS, RS, tau=NULL){
  if(is.null(tau)){
    out = sort(unique(c(LE, RE, LS, RS)))
  }else{
    out = sort(unique(c(LE, RE, LS, RS, tau)))
  }
}


#' Map Breakpoints to Their Ranks
#'
#' Map each endpoint to its corresponding position in a sorted
#' vector of breakpoints.
#'
#' @param LE A vector of left endpoints.
#' @param RE A vector of right endpoints.
#' @param LS A vector of left endpoints for the second set.
#' @param RS A vector of right endpoints for the second set.
#' @param breaks A sorted vector of unique breakpoints against
#'   which the endpoints are matched.
#'
#' @return A named list with four components, \code{LE}, \code{RE},
#'   \code{LS}, and \code{RS}, containing the integer indices
#'   of the corresponding endpoints in \code{breaks}.
#'   Endpoints not found in \code{breaks} are returned as \code{NA}.
#'
#' @keywords internal
rankbreaks <- function(LE, RE, LS, RS, breaks){
  list(LE=match(LE, breaks),
       RE=match(RE, breaks),
       LS=match(LS, breaks),
       RS=match(RS, breaks))
}


#' Fit a Conditional Model
#'
#' Fit the conditional model using variational Bayes (VB) or the
#' expectation-maximization (EM) algorithm.
#'
#' The breakpoints are constructed from the unique values of
#' \code{LE}, \code{RE}, \code{LS}, and \code{RS}.
#' The endpoints are then mapped to their ranks in the breakpoint
#' vector and passed to the underlying C++ implementation.
#'
#' @param LE A numeric vector of left endpoints for the first set
#'   of observations.
#' @param RE A numeric vector of right endpoints for the first set
#'   of observations.
#' @param LS A numeric vector of left endpoints for the second set
#'   of observations.
#' @param RS A numeric vector of right endpoints for the second set
#'   of observations.
#' @param maxit An integer specifying the maximum number of
#'   iterations. Defaults to \code{1000L}.
#' @param tol A positive numeric value specifying the convergence
#'   tolerance for the iterative algorithm. Defaults to \code{1e-5}.
#' @param method A character string specifying the estimation method.
#'   Either \code{"VB"} for variational Bayes or \code{"EM"} for
#'   the expectation-maximization algorithm. Defaults to \code{"VB"}.
#'
#' @return A list containing the results returned by the underlying
#'   estimation procedure, with the following modifications:
#'   \describe{
#'     \item{logprob}{The log-probability vector with its first
#'       element removed.}
#'     \item{breaks}{The sorted vector of unique breakpoints
#'       constructed from the input endpoints.}
#'   }
#'   Other components depend on the underlying estimation procedure.
#'
#' @seealso \code{\link{conditional_dic_gibbs}},
#'   \code{\link{jointfit_dic}}, \code{\link{joint_dic_gibbs}}
#'
#' @export
conditionalfit_dic <- function(LE, RE, LS, RS,
                               maxit=1000L, tol=1e-5, method="VB"){
  breaks <- setbreaks(LE, RE, LS, RS)
  RANK <- rankbreaks(LE, RE, LS, RS, breaks)
  #the index is 0-origin in cpp
  if(method == "VB"){
    res <- with(RANK, conditional_DIC_vb(LE-1L, RE-1L, LS-1L, RS-1L, 
                                   m = length(breaks), maxit = maxit, tol = tol))
  }else if(method == "EM"){
     res <- with(RANK, conditional_DIC_em(LE-1L, RE-1L, LS-1L, RS-1L, 
                                    m = length(breaks), maxit = maxit, tol = tol))    
  }else{
    warning("this method is not implemented")
    #break
  }
  res$logprob <- res$logprob[-1]
  res$breaks <- breaks
  return(res)
}


#' Fit a Conditional Model Using Gibbs Sampling
#'
#' Fit the conditional model using Gibbs sampling and compute
#' the deviance information criterion (DIC).
#'
#' The breakpoints are constructed from the unique values of
#' \code{LE}, \code{RE}, \code{LS}, and \code{RS}.
#' The endpoints are mapped to their ranks in the breakpoint
#' vector and passed to the underlying C++ implementation.
#'
#' @param LE A numeric vector of left endpoints for the first set
#'   of observations.
#' @param RE A numeric vector of right endpoints for the first set
#'   of observations.
#' @param LS A numeric vector of left endpoints for the second set
#'   of observations.
#' @param RS A numeric vector of right endpoints for the second set
#'   of observations.
#' @param maxit An integer specifying the number of Gibbs sampling
#'   iterations. Defaults to \code{1000L}.
#'
#' @return A list containing the results returned by the underlying
#'   Gibbs sampling procedure, with an additional component:
#'   \describe{
#'     \item{breaks}{The sorted vector of unique breakpoints
#'       constructed from the input endpoints.}
#'   }
#'   Other components depend on the underlying estimation procedure.
#'
#' @seealso \code{\link{conditionalfit_dic}},
#'   \code{\link{jointfit_dic}}, \code{\link{joint_dic_gibbs}}
#'
#' @export
conditional_dic_gibbs <- function(LE, RE, LS, RS, iter=1000L){
  breaks <- setbreaks(LE, RE, LS, RS)
  RANK <- rankbreaks(LE, RE, LS, RS, breaks)
  res <- with(RANK, conditional_DIC_gibbs(LE-1L, RE-1L, LS-1L, RS-1L, 
                                          m = length(breaks),
                                          iter = iter))
  res$breaks <- breaks
  return(res)
}


#' Fit a Joint Model and Compute the Deviance Information Criterion
#'
#' Fit the joint model using variational Bayes (VB) or the
#' expectation-maximization (EM) algorithm, and compute the
#' deviance information criterion (DIC).
#'
#' The breakpoints are constructed from the unique values of
#' \code{LE}, \code{RE}, \code{LS}, and \code{RS}, together with
#' the optional additional breakpoints in \code{tau}.
#' The endpoints are mapped to their ranks in the breakpoint
#' vector and passed to the underlying C++ implementation.
#'
#' @param LE A numeric vector of left endpoints for the first set
#'   of observations.
#' @param RE A numeric vector of right endpoints for the first set
#'   of observations.
#' @param LS A numeric vector of left endpoints for the second set
#'   of observations.
#' @param RS A numeric vector of right endpoints for the second set
#'   of observations.
#' @param tau An optional numeric vector of additional breakpoints.
#'   Defaults to \code{NULL}, in which case no additional
#'   breakpoints are included.
#' @param maxit An integer specifying the maximum number of
#'   iterations. Defaults to \code{1000L}.
#' @param tol A positive numeric value specifying the convergence
#'   tolerance for the iterative algorithm. Defaults to \code{1e-5}.
#' @param method A character string specifying the estimation method.
#'   Either \code{"VB"} for variational Bayes or \code{"EM"} for
#'   the expectation-maximization algorithm. Defaults to \code{"VB"}.
#'
#' @return A list containing the results returned by the underlying
#'   estimation procedure, with the following modifications:
#'   \describe{
#'     \item{logprob}{The log-probability vector with its first
#'       element removed.}
#'     \item{breaks}{The sorted vector of unique breakpoints
#'       constructed from the input endpoints and, if supplied,
#'       \code{tau}.}
#'   }
#'   Other components depend on the underlying estimation procedure.
#'
#' @seealso \code{\link{joint_dic_gibbs}},
#'   \code{\link{conditionalfit_dic}},
#'   \code{\link{conditional_dic_gibbs}}
#'
#' @export
jointfit_dic <- function(LE, RE, LS, RS, tau=NULL, maxit=1000L, tol=1e-5, method="VB"){
  breaks <- setbreaks(LE, RE, LS, RS, tau=tau)
  RANK <- rankbreaks(LE, RE, LS, RS, breaks)
  #the index is 0-origin in cpp
  if(method == "VB"){
    res <- with(RANK, joint_DIC_vb(LE-1L, RE-1L, LS-1L, RS-1L, 
                                   m = length(breaks), maxit = maxit, tol = tol))
  }else if(method == "EM"){
    res <- with(RANK, joint_DIC_em(LE-1L, RE-1L, LS-1L, RS-1L, 
                                   m = length(breaks), maxit = maxit, tol = tol))    
  }else{
    warning("this method is not implemented")
    #break
  }
  res$logprob <- res$logprob[-1]
  res$breaks <- breaks
  return(res)
}


#' Fit a Joint Model Using Gibbs Sampling
#'
#' Fit the joint model using Gibbs sampling and compute
#' the deviance information criterion (DIC).
#'
#' The breakpoints are constructed from the unique values of
#' \code{LE}, \code{RE}, \code{LS}, and \code{RS}, together with
#' the optional additional breakpoints in \code{tau}.
#' The endpoints are mapped to their ranks in the breakpoint
#' vector and passed to the underlying C++ implementation.
#'
#' @param LE A numeric vector of left endpoints for the first set
#'   of observations.
#' @param RE A numeric vector of right endpoints for the first set
#'   of observations.
#' @param LS A numeric vector of left endpoints for the second set
#'   of observations.
#' @param RS A numeric vector of right endpoints for the second set
#'   of observations.
#' @param tau An optional numeric vector of additional breakpoints.
#'   Defaults to \code{NULL}, in which case no additional
#'   breakpoints are included.
#' @param iter An integer specifying the number of Gibbs sampling
#'   iterations. Defaults to \code{1000L}.
#'
#' @return A list containing the results returned by the underlying
#'   Gibbs sampling procedure, with the following modifications:
#'   \describe{
#'     \item{logprob}{The log-probability vector with its first
#'       element removed.}
#'     \item{breaks}{The sorted vector of unique breakpoints
#'       constructed from the input endpoints and, if supplied,
#'       \code{tau}.}
#'   }
#'   Other components depend on the underlying estimation procedure.
#'
#' @seealso \code{\link{jointfit_dic}},
#'   \code{\link{conditionalfit_dic}},
#'   \code{\link{conditional_dic_gibbs}}
#'
#' @export
joint_dic_gibbs <- function(LE, RE, LS, RS, tau=NULL,iter=1000L){
  breaks <- setbreaks(LE, RE, LS, RS, tau=tau)
  RANK <- rankbreaks(LE, RE, LS, RS, breaks)
  #the index is 0-origin in cpp
  res <- with(RANK, joint_DIC_gibbs(LE-1L, RE-1L, LS-1L, RS-1L, 
                                    m = length(breaks),
                                    iter = iter))
  res$logprob <- res$logprob[-1]
  res$breaks <- breaks
  return(res)
}

