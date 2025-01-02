colMarginal <- function(X){
  cumsum(colSums(X))  
}

p2ccdf <- function(p){
  rev(cumsum(rev(p)))
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


# revmat <- function(rand){
#   apply(rand, 2, rev)
# }