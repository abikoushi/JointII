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
set.seed(12345)
at = weibull_process(tau, shape = sigma, maxit = 10000)
y <- rweibull(length(at), trueparam$shape, scale = trueparam$scale)
mean(y)
simDat0 <- simDIC(y, at = at, WS = 1, WE = 1)
simDat <- simDat0

system.time({
  res <- conditional_dic_gibbs(LE = simDat$LE, RE = simDat$RE, LS = simDat$LS,　RS =  simDat$RS,
                               maxit=100L)
})

plot(res$logprob, type = "l")
Fhat = rowMeans(1-apply(res$incubation[,51:100],2,cumsum))
plot(res$breaks, Fhat, type="s", col="royalblue", xlim = c(0,40))
curve(pweibull(x, shape = trueparam$shape, scale = trueparam$scale, lower.tail = FALSE), add=TRUE, lty=2)

####

system.time({
  est1 <- conditionalfit_dic(LE = simDat$LE, RE = simDat$RE, LS = simDat$LS,　RS =  simDat$RS, 
                             maxit=100L, tol=1e-8, method = "EM")
})

system.time({
  est2 <- conditionalfit_dic(LE = simDat$LE, RE = simDat$RE, LS = simDat$LS,　RS =  simDat$RS,
                             maxit=100L, tol=1e-8, method = "VB")
})

est2$logprob
#wch = which(est2$alpha > .Machine$double.eps)

plot(est1$logprob[-1], type = "l")
plot(est2$logprob[-1], type = "l")

plot(est1$breaks, prob2ccdf(est1$incubation), type="s", col="royalblue", xlim = c(0,40))
lines(est2$breaks, prob2ccdf(est2$incubation), type="s", col="forestgreen", xlim = c(0,40))
lines(res$breaks, Fhat, col="orange", type="s")
curve(pweibull(x, shape = trueparam$shape, scale = trueparam$scale, lower.tail = FALSE), add=TRUE, lty=2)

