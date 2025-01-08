library(JointII)

trueparam <- list(shape=2, scale = 1)
integ <- function(x){
  pweibull(tau-x, trueparam[[1]], trueparam[[2]], lower.tail=FALSE) * lambda(x)  
}
lambda <- function(x){
  sigma*(x^(sigma-1))
}
Lambda <- function(x){
  x^sigma
}

sigma = 0.8
res_int <- integrate(integ, lower = 0, upper = tau)
Ns <- integer(1e+4L)
for(i in 1L:1e+4L){
  at <- weibull_process(window = tau, shape = sigma, maxit = 100000)
  y <- rweibull(length(at), trueparam[[1]], trueparam[[2]])
  Ns[i] <- sum(at+y>tau)
}

round(mean(Ns),2)
round(res_int$value,2)
tab <- table(Ns)
y <- as.integer(names(tab))
plot(y, as.numeric(tab)/10000, type = "b")
points(y, dpois(y, lambda = res_int$value), pch=2, type = "b")



tau = 20
sigma = 1.2
res <- vector("list", 100)
for(i in 1:100){
  res[[i]] <- weibull_process(window = tau, shape = sigma, maxit = 100000)  
}

curve(x^sigma, 0, tau, ylim=c(0,60))
for(i in 1:100){
  lines(res[[i]], 1:length(res[[i]]), type = "s", col=rgb(0,0,0,0.1))
}


integrate(integ, lower = 0, upper = tau)
trueparam$scale*gamma(1+1/trueparam$shape)

