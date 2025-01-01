p2ccdf <- function(p){
  rev(cumsum(rev(p)))
}

colMarginal <- function(X){
  cumsum(colSums(X))  
}

IIdf <- function(est){
  data.frame(
    time = est$breaks,
    ccdf = p2ccdf(est$incubation),
    intensity = cumsum(est$infection)*(sum(est$event)+sum(est$truncated_event)),
    cuminf = colMarginal(est$event),
    cuminf_unobs = colMarginal(est$truncated_event)
  )
}

rmlast <- function(x)x[-length(x)]

setbreaks <- function(LE, RE, LS, RS){
  sort(unique(c(LE, RE, LS, RS)))
}

rankbreaks <- function(LE, RE, LS, RS, breaks){
  list(LE=match(LE, breaks),
       RE=match(RE, breaks),
       LS=match(LS, breaks),
       RS=match(RS, breaks))
}

jointecdf_dic_em <- function(LE, RE, LS, RS, maxit=1000L, tol=1e-5){
  breaks <- setbreaks(LE, RE, LS, RS)
  RANK <- rankbreaks(LE, RE, LS, RS, breaks)
  #the index is 0-origin in cpp
  res <- with(RANK, joint_DIC_em(LE-1L, RE-1L, LS-1L, RS-1L, 
                      m = length(breaks), maxit = maxit, tol = tol))
  res$logprob <- res$logprob[-1]
  res$breaks <- breaks
  return(res)
}

jointecdf_dic_vb <- function(LE, RE, LS, RS, maxit=1000L, tol=1e-5){
  breaks <- setbreaks(LE, RE, LS, RS)
  RANK <- rankbreaks(LE, RE, LS, RS, breaks)
  #the index is 0-origin in cpp
  res <- with(RANK, joint_DIC_vb(LE-1L, RE-1L, LS-1L, RS-1L, 
                                 m = length(breaks), maxit = maxit, tol = tol))
  res$logprob <- res$logprob[-1]
  res$breaks <- breaks
  return(res)
}