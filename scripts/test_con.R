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
#library(RColorBrewer)

juliapath <- scan("path.txt", what = character())
julia <- julia_setup(JULIA_HOME = juliapath[1])
julia_source(juliapath[4])
julia_source(juliapath[5])


Lambda <- function(x, sigma){
  x^sigma
}

trueparam <- list(shape=1.5, scale=5)
sigma <- 1
tau <- 100L
set.seed(111)
at = weibull_process(tau, shape = sigma, maxit = 10000)
y <- rweibull(length(at), trueparam$shape, scale = trueparam$scale)
mean(y)
simDat0 <- simDIC(y, at = at, WS = 1, WE = 1)
simDat <- simDat0

system.time({
  res <- conditional_dic_gibbs(LE = simDat$LE, RE = simDat$RE, LS = simDat$LS, RS =  simDat$RS,
                               maxit=2000L)
})

res2 <- julia_call("empest.gibbssampler", simDat$LE, simDat$RE, simDat$LS, simDat$RS, 2000L)

plot(res2$logprob, type = "l")
lines(res$logprob, type = "l", col="blue")
burnin <- 1:100
Shat = rowMeans(1-apply(res$incubation[,-burnin],2,cumsum))
Shat2 = colMeans(res2$ccdf[-burnin,])

plot(res$breaks, Shat, type="s", col="royalblue", xlim = c(0,40), ylim = c(0,1))
lines(res$breaks, Shat2, type="s", col="orange")
curve(pweibull(x, shape = trueparam$shape, scale = trueparam$scale, lower.tail = FALSE), add=TRUE, lty=2)

####

system.time({
  est_vb <- conditionalfit_dic(LE = simDat$LE, RE = simDat$RE, LS = simDat$LS,　RS =  simDat$RS, 
                             maxit=200L, tol=1e-8, method = "VB")
})

system.time({
  est_vb2 <- julia_call("contool.vb", simDat$LE, simDat$RE, simDat$LS, simDat$RS, 200L)
})

plot(est_vb$logprob, type = "l")
lines(est_vb2$logprob, type = "l", col="red")

plot(est_vb$breaks, prob2ccdf(est_vb$incubation), type="s", col="royalblue", xlim = c(0,40))
lines(est_vb2$value, prob2ccdf(est_vb2$alpha/sum(est_vb2$alpha)), type="s", col="orange", xlim = c(0,40))
curve(pweibull(x, shape = trueparam$shape, scale = trueparam$scale, lower.tail = FALSE), add=TRUE, lty=2)

system.time({
  est_em <- conditionalfit_dic(LE = simDat$LE, RE = simDat$RE, LS = simDat$LS,　RS =  simDat$RS,
                             maxit=200L, tol=1e-8, method = "EM")
})

system.time({
  est_em2 <- julia_call("contool.EM", simDat$LE, simDat$RE, simDat$LS, simDat$RS,200L)
})

plot(est_em$logprob[-1], type = "l")
lines(est_em2$logprob[-1], type = "l", col="red")

plot(est_vb$breaks, prob2ccdf(est_em$incubation), type="s", col="royalblue", xlim = c(0,40))
lines(est_vb2$value, prob2ccdf(est_em2$alpha/sum(est_em2$alpha)), type="s", col="orange", xlim = c(0,40))
curve(pweibull(x, shape = trueparam$shape, scale = trueparam$scale, lower.tail = FALSE), add=TRUE, lty=2)

