rmlast <- function(x)x[-length(x)]

setbreaks <- function(LE, RE, LS, RS, tau=NULL){
  if(is.null(tau)){
    out = sort(unique(c(LE, RE, LS, RS)))
  }else{
    out = sort(unique(c(LE, RE, LS, RS, tau)))
  }
}

rankbreaks <- function(LE, RE, LS, RS, breaks){
  list(LE=match(LE, breaks),
       RE=match(RE, breaks),
       LS=match(LS, breaks),
       RS=match(RS, breaks))
}

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

conditional_dic_gibbs <- function(LE, RE, LS, RS, maxit=1000L){
  breaks <- setbreaks(LE, RE, LS, RS)
  RANK <- rankbreaks(LE, RE, LS, RS, breaks)
  res <- with(RANK, conditional_DIC_gibbs(LE-1L, RE-1L, LS-1L, RS-1L, 
                                          m = length(breaks),
                                          iter = maxit))
  res$breaks <- breaks
  return(res)
}

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

normalize_cdf <- function(rand){
  sweep(rand,2,rand[nrow(rand),],"/")
}

rand_cdf <- function(alpha, np, normalize=FALSE){
  out = rand_gamma(alpha, np=np)
  out = apply(out, 1, cumsum)
  if(normalize){
    out = normalize_cdf(out)
  } 
  return(out)
}

prob2ccdf <- function(h){
  rev(cumsum(rev(h)))
}
