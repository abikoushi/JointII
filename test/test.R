library(JointII)
library(Matrix)
library(bench)
library(ggplot2)
library(pammtools)
#library(RColorBrewer)
#weibull_process(100, shape = 1.5, maxit = 1000)
y <- rweibull(50, 2, 5)
simDat <- simDIC(y)
head(simDat)
bmrec <- bench::mark({
  est <- jointecdf_dic_vb(simDat$LE, simDat$RE, simDat$LS, simDat$RS,
                          maxit=1000L, tol=1e-5)
}, min_iterations = 1, max_iterations = 1)

plot(est$logprob, type="l")

revmat <- function(rand){
  apply(rand, 2, rev)
}

normalize_cdf <- function(rand){
  sweep(rand,2,rand[nrow(rand),],"/")
}
rand_cum
rand_cdf <- function(alpha, np=2000L){
  wch = which(alpha>0)
  out = matrix(0, nrow=np, ncol=length(alpha))
  for(i in wch){
    out[,i] = rgamma(np, alpha[i], 1)    
  }
  apply(out, 1, cumsum)
}



k01 = rgb(0,0,0,0.01)

plot(est$breaks,cumsum(est$alpha_l), type="s", col="royalblue", ylim=c(0,100),xlim=c(0,100))
abline(0,1,col="orangered")
rand = rand_cdf(est$alpha_l)
matlines(est$breaks, rand, col=k01, lty=1)

curve(50*pweibull(x,2,5), type="l", col="orangered", xlim = range(est$breaks),ylim=c(0,70))
rand = rand_cdf(est$alpha_h, np = 1000)
matlines(est$breaks, rand, col=k01, lty=1)

dim(rand)
rand = rand_cdf(est$alpha_h, np = 1000)
curve(pweibull(x,2,5,lower.tail = TRUE), type="l", col="orangered", xlim = range(est$breaks))
matlines(est$breaks, normalize_cdf(rand), type="s", lty = 1,col=k01)

