
simDIC <- function(y, at=NULL, WS=NULL, WE=NULL){
  N <- length(y)
  if(is.null(at)){
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
  r_s <- runif(N)
  LS <- S - r_s*WS
  RS <- S + (1-r_s)*WS
  LE <- at - r_e*WE
  RE <- at + (1-r_e)*WE
  return(data.frame(LE=pmax(LE,0), RE=RE,
                    LS=pmax(LS,LE,0), RS=RS,
                    E=at, S=S))
}

weibull_process <- function(window, shape, maxit){
  out = numeric(maxit)
  at = 0.0
  z = runif(1)
  at = (-log(z))^(1/shape)
  out[1] = at
  for(i in 1:maxit){
    z = runif(1) #next infection
    at =  ((-log(z))^(1/shape) + at^shape)^(1/shape)
    out[i] = at
    if(at >= window){
      out = out[1:i]
      break
    }
  }
  return(out)
}
