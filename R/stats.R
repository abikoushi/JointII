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

simulated_pvalue <-function(rand, est, tp, brks, tvalue){
  CI <- matrix(0, nrow = length(brks), ncol = 3)
  colnames(CI) <- c("lower", "upper", "estimates")
  
  xv = numeric(length(brks))
  pv <- numeric(length = length(brks))
  for(i in 1:length(brks)){
    ind = which.max(est$breaks > brks[i])-1
    ps <- seq(0, 1, length.out=ncol(rand))
    rs <- sort(rand[ind,])
    lower <- ps[which.min(rs < tvalue[i])]
    upper <- rev(ps)[which.min(tvalue[i]<rev(rs))]
    pv[i] <- min(lower, upper)
    
    xv[i] <- est$breaks[ind]
    ci <- quantile(rand[ind,], probs = c(pv[i], 1-pv[i], 0.5))
    CI[i,] <- ci
  }
  data.frame(x=brks, pvalue=pv, CI)
}

# revmat <- function(rand){
#   apply(rand, 2, rev)
# }