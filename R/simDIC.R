#' Simulate Interval-Censored Data for DIC Analysis
#'
#' Generate a simulated dataset consisting of interval-censored
#' observations for two event times.
#'
#' Event times are generated from the supplied observation times
#' and interval widths. If the observation times or interval widths
#' are not supplied, they are generated from exponential
#' distributions.
#'
#' @param y A numeric vector of non-negative durations between
#'   the first and second event times.
#' @param at An optional numeric vector of observation times for
#'   the first event. If \code{NULL}, the observation times are
#'   generated as cumulative sums of independent exponential
#'   random variables, corresponding to the arrival times of a
#'   homogeneous Poisson process with unit rate.
#'   Defaults to \code{NULL}.
#' @param WS An optional numeric vector of interval widths for
#'   the second event. If \code{NULL}, the widths are generated
#'   independently from an exponential distribution with unit rate.
#'   Defaults to \code{NULL}.
#' @param WE An optional numeric vector of interval widths for
#'   the first event. If \code{NULL}, the widths are generated
#'   independently from an exponential distribution with unit rate.
#'   Defaults to \code{NULL}.
#'
#' @return A data frame with \code{length(y)} rows and six columns:
#'   \describe{
#'     \item{LE}{The lower endpoint of the observation interval
#'       for the first event, truncated at zero.}
#'     \item{RE}{The upper endpoint of the observation interval
#'       for the first event.}
#'     \item{LS}{The lower endpoint of the observation interval
#'       for the second event, truncated at zero and constrained
#'       to be no smaller than \code{LE}.}
#'     \item{RS}{The upper endpoint of the observation interval
#'       for the second event.}
#'     \item{E}{The generated time of the first event.}
#'     \item{S}{The generated time of the second event,
#'       obtained as \code{at + y}.}
#'   }
#'
#' @details
#' For each observation, the interval endpoints are generated
#' by randomly positioning the event time within an interval
#' of width \code{WE} or \code{WS}, respectively.
#' The lower endpoints are truncated at zero, and the lower
#' endpoint of the second event is additionally constrained
#' to be no smaller than the lower endpoint of the first event.
#'
#' @seealso \code{\link{weibull_process}}
#'
#' @export
simDIC <- function(y, at=NULL, WS=NULL, WE=NULL){
  N <- length(y)
  if(is.null(at)){
    #homogineous Poisson process
    at <- cumsum(rexp(N))
  }
  if(is.null(WS)){
    WS <- rexp(N)
  }
  if(is.null(WE)){
    WE <- rexp(N)
  }
  S <- at + y
  r_e <- runif(N)
  LE <- at - r_e*WE
  RE <- at + (1-r_e)*WE
  r_s <- runif(N)
  LS <- S - r_s*WS
  RS <- S + (1-r_s)*WS
  return(data.frame(LE=pmax(LE,0), RE=RE,
                    LS=pmax(LS,LE,0), RS=RS,
                    E=at, S=S))
}


#' Generate Arrival Times from a Weibull Process
#'
#' Generate event arrival times from a renewal process whose
#' interarrival times follow a Weibull distribution with
#' the specified shape parameter.
#'
#' The process is generated sequentially until the next event
#' time reaches or exceeds the specified observation window.
#'
#' @param window A positive numeric value specifying the end
#'   of the observation window.
#' @param shape A positive numeric value specifying the shape
#'   parameter of the Weibull distribution.
#' @param maxit An integer specifying the maximum number of
#'   event times to allocate and generate. Defaults to the
#'   value supplied by the user.
#'
#' @return A numeric vector of event arrival times occurring
#'   before \code{window}. The returned times are in increasing
#'   order.
#'
#' @details
#' The first event time is generated as
#' \deqn{T_1 = (-\log U_1)^{1/\mathrm{shape}},}
#' where \eqn{U_1} is uniformly distributed on \eqn{(0,1)}.
#'
#' Subsequent event times are generated recursively as
#' \deqn{T_i = \left(-\log U_i + T_{i-1}^{\mathrm{shape}}
#' \right)^{1/\mathrm{shape}},}
#' where the \eqn{U_i} are independent uniform random variables.
#'
#' This construction corresponds to a non-homogeneous Poisson
#' process with cumulative intensity
#' \eqn{\Lambda(t)=t^{\mathrm{shape}}}.
#'
#' @seealso \code{\link{simDIC}}
#'
#' @export
weibull_process <- function(window, shape, maxit){
  out = numeric(maxit)
  at = 0.0
  z = runif(1)
  at = (-log(z))^(1/shape)
  out[1] = at
  for(i in 1:maxit){
    z = runif(1) #next infection
    at =  (-log(z) + at^shape)^(1/shape)
    #at =  ((-log(z))^(1/shape) + at^shape)^(1/shape)
    if(at >= window){
      out = out[1:(i-1)]
      break
    }
    out[i] = at
  }
  return(out)
}
