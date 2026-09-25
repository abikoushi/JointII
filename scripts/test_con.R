library(JointII)
library(Matrix)
library(bench)
library(ggplot2)
library(pammtools)
library(reshape2)
library(gridExtra)
library(ggplot2)
library(dplyr)
#library(RColorBrewer)

juliapath <- scan("path.txt", what = character())
julia <- julia_setup(JULIA_HOME = juliapath[1])
julia_source(juliapath[4])

logsumexp <- function(x){
  mx <- max(x)
  mx + log(sum(exp(x-mx)))
}

# logsumexp(1:10)
# log(sum(exp(1:10)))

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
  res <- conditional_dic_gibbs(LE = simDat$LE, RE = simDat$RE, LS = simDat$LS,　RS =  simDat$RS,
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
  est1 <- conditionalfit_dic(LE = simDat$LE, RE = simDat$RE, LS = simDat$LS,　RS =  simDat$RS, 
                             maxit=200L, tol=1e-8, method = "EM")
})

system.time({
  est2 <- conditionalfit_dic(LE = simDat$LE, RE = simDat$RE, LS = simDat$LS,　RS =  simDat$RS,
                             maxit=200L, tol=1e-8, method = "VB")
})


plot(est1$logprob[-1], type = "l")
lines(est2$logprob[-1], type = "l")

plot(est1$breaks, prob2ccdf(est1$incubation), type="s", col="royalblue", xlim = c(0,40))
lines(est2$breaks, prob2ccdf(est2$incubation), type="s", col="forestgreen", xlim = c(0,40))
lines(res$breaks, Fhat, col="orange", type="s")
curve(pweibull(x, shape = trueparam$shape, scale = trueparam$scale, lower.tail = FALSE), add=TRUE, lty=2)

