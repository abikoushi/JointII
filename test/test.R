library(JointII)
library(Matrix)
library(bench)
library(ggplot2)
library(pammtools)
#library(RColorBrewer)

#weibull_process(100, shape = 1.5, maxit = 1000)
y <- rweibull(50, 2, 5)
simDat <- simDIC(y)

bmrec <- bench::mark({
  est <- jointfit_dic(simDat$LE, simDat$RE, simDat$LS, simDat$RS,
                          maxit=1000L, tol=1e-5)
}, min_iterations = 1, max_iterations = 1)

plot(est$logprob, type="l")

k01 = rgb(0,0,0,0.01)

plot(est$breaks,cumsum(est$beta), type="s", col="royalblue", ylim=c(0,100),xlim=c(0,100))
abline(0,1,col="orangered")
rand = rand_cdf(est$beta)
matlines(est$breaks, rand, col=k01, lty=1)

#png("ccdf.png")
par(mfrow=c(2,1))
curve(50*pweibull(x,2,5), type="l", col="orangered", xlim = range(est$breaks),ylim=c(0,80))
rand = rand_cdf(est$alpha, np = 1000)
matlines(est$breaks, rand, col=k01, lty=1)
curve(pweibull(x,2,5), type="l", col="orangered", xlim = range(est$breaks))
matlines(est$breaks, normalize_cdf(rand), type="s", lty = 1,col=k01)
#dev.off()
