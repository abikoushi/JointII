library(JointII)
library(Matrix)

Lambda <- function(x, sigma){
  x^sigma
}


#trueparam <- list(shape=1.5, scale=5)
trueparam <- list(meanlog=log(10) - 0.5, sdlog=1)
sigma <- 1
tau <- 100L
prob <- c(0.5, 0.8, 0.9, 0.95)
qs <- qlnorm(prob, trueparam[[1]], trueparam[[2]])

pv_s <- matrix(0,100,4)
for(it in 1:100){
  set.seed(it)
  at = weibull_process(tau, shape = sigma, maxit = 10000)
  y <- rlnorm(length(at), trueparam[[1]], trueparam[[2]])
  
  simDat0 <- simDIC(y, at = at, WS = 2, WE = 2)
  simDat <- simDat0[simDat0$RS < tau,]
  
  system.time({
    res <- joint_dic_gibbs(LE = simDat$LE, RE = simDat$RE, LS = simDat$LS,　RS =  simDat$RS, tau = tau,
                           iter=4000L)
  })
  
  burnin <- 1:2000
  #curve(plnorm(x, trueparam[[1]], trueparam[[2]], lower.tail = FALSE), 0, 50, ylab = "ccdf")
  #lines(res$breaks, Shat)
  
  # Stilde <- apply(res$incubation[,-burnin],2,prob2ccdf)
  # 
  # matlines(res$breaks, Stilde, lwd=0.1, col=rgb(0,0,0,0.1))
  # Shat <- rowMeans(apply(res$incubation[,-burnin],2,prob2ccdf))
  
  # curve(plnorm(x, trueparam[[1]], trueparam[[2]], lower.tail = FALSE), 0, 50, ylab = "ccdf")
  # for(i in 1:2000){
  # ind <- which(res$incubation[,-burnin][,i] > 0.0)
  # afun <- approx(x = res$breaks[ind], y = prob2ccdf(res$incubation[ind,-burnin][,i]), n = 100)
  # lines(afun$x, afun$y, type = "l", lwd=0.1, col=rgb(0,0,0,0.1))
  # }
  # 
  
  ptilde <- matrix(0, 2000, 4)
  for(i in 1:2000){
    ind <- which(res$incubation[,-burnin][,i] > 0.0)
    afun <- approx(x = res$breaks[ind], y = prob2ccdf(res$incubation[ind,-burnin][,i]), xout = qs)
    ptilde[i,] <- ifelse(is.na(afun$y),0,afun$y)
  }
  
  for(i in 1:4){
    lower <- ecdf(ptilde[,i])(1-prob[i])
    pv_s[it, i] <- 2*min(lower, 1-lower)  
  }  
}


