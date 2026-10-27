library(JointII)
library(Matrix)
library(bench)
library(ggplot2)
library(pammtools)
library(reshape2)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(JuliaCall)
juliapath <- scan("path.txt", what = character())
julia <- julia_setup(JULIA_HOME = juliapath[1])
julia_source(juliapath[2])
julia_source(juliapath[3])


Lambda <- function(x, sigma){
  x^sigma
}


trueparam <- list(shape=1.5, scale=10)
sigma <- 1.0
tau <- 100L
set.seed(111)
curve(Lambda(x,sigma),0,100)
at = weibull_process(tau, shape = sigma, maxit = 10000)
lines(at[-1], 1:length(at[-1]), type = "s", col="grey")

y <- rweibull(length(at), trueparam$shape, scale = trueparam$scale)

WE=1
WS=1

simDat0 <- simDIC(y, at = at, WS = 1, WE = 1)
simDat <- simDat0[simDat0$RS < tau,]
system.time({
  res <- joint_dic_gibbs(LE = simDat$LE, RE = simDat$RE, LS = simDat$LS, RS =  simDat$RS, tau = tau,
                         iter=2000L)
})

res2 <- julia_call("jointest.gibbssampler", simDat$LE, simDat$RE, simDat$LS, simDat$RS, tau, 2000L)

burnin <- 1:100
plot(res$logprob, type = "l")
lines(res2$logprob, type = "l", col="cornflowerblue")

Shat = rowMeans(1-apply(res$incubation[,-burnin],2,cumsum))
Shat2 = colMeans(res2$ccdf[-burnin,])



plot(res$breaks, Shat, type="s", col="royalblue", xlim = c(0,100))
lines(res2$value, Shat2, type="s", col="orange", xlim = c(0,100))
curve(pweibull(x, shape = trueparam$shape, scale = trueparam$scale, lower.tail = FALSE), add=TRUE, lty=2)


Lamhat = rowMeans(res$infection[,51:100])
Lamhat2 = colMeans(res2$intensity[-burnin,])
plot(res$breaks, Lamhat, type="s", col="royalblue", xlim = c(0,100))
lines(res2$value, Lamhat2, type="s", col="orange", xlim = c(0,100))
abline(0,1)


system.time({
  est_vb <- jointfit_dic(LE = simDat$LE, RE = simDat$RE, LS = simDat$LS,　RS =  simDat$RS, tau = tau,
                      maxit=200L, tol=1e-8, method = "VB")
})
system.time({
  est_vb2 <- julia_call("jointtool.vb", simDat$LE, simDat$RE, simDat$LS, simDat$RS, tau, 200L)
})

plot(est_vb$logprob, type = "l")
lines(est_vb2$logprob, type = "l")
#plot(est_vb$breaks, cumsum(est_vb$incubation), type="s")
plot(est_vb$breaks, cumsum(est_vb$alpha/sum(est_vb$alpha)), type="s")
lines(est_vb2$value, cumsum(est_vb2$alpha/sum(est_vb2$alpha)), type="s", col="red")

plot(est_vb$breaks, cumsum(est_vb$alpha/sum(est_vb$alpha)), type="s")
lines(est_vb2$value, cumsum(est_vb2$alpha/sum(est_vb2$alpha)), type="s", col="red")

plot(est_vb$breaks, cumsum(est_vb$beta), type="s")
lines(est_vb2$value, cumsum(est_vb2$beta), type="s", col="red")


system.time({
  est_em <- jointfit_dic(LE = simDat$LE, RE = simDat$RE, LS = simDat$LS,　RS =  simDat$RS, tau = tau,
                      maxit=1000L, tol=1e-8, method = "EM")
})
system.time({
  est_em2 <- julia_call("jointtool.EM", simDat$LE, simDat$RE, simDat$LS, simDat$RS, tau)  
})

plot(est_em$breaks, cumsum(est_em$incubation), type="s")
lines(est_em2$value, cumsum(est_em2$alpha/sum(est_em2$alpha)), type="s", col="red")

plot(est_em$breaks, cumsum(est_em$beta), type="s")
lines(est_em2$value, cumsum(est_em2$beta), type="s", col="red")

