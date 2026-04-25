library(JointII)
library(Matrix)
library(bench)
library(ggplot2)
library(pammtools)
library(reshape2)
library(gridExtra)
library(ggplot2)
library(dplyr)

# x <- rgamma(10,1)
# x <- x/sum(x)
# X <- sparseMatrix(i = 1:10, j = 1:10, x = x)
# res <- JointII:::rcate(10000,X)
# 
# plot(diag(X),diag(res)/sum(diag(res)))
# abline(0,1)

logsumexp <- function(x){
  mx <- max(x)
  mx + log(sum(exp(x-mx)))
}

weibull_process = function(tau, shape, maxit){
  w = numeric(maxit)
  z = runif(1)
  at =  (-log(z))^(1/shape) 
  if(maxit >= 2){
    for(i in 1:maxit){
      z = runif(1)
      at =  (-log(z) + at^shape)^(1/shape)
      w[i] = at
      if(at>tau){
        break
      }
    }    
  }
  return(w[1:i])
}

Lambda <- function(x, sigma){
  x^sigma
}
simDIC <- function(y,at,WS,WE){
  S <- at+y
  ue = runif(length(at))
  us = runif(length(at))
  LE = pmax(0, at - WE * (1 - ue))
  LS = pmax(LE, S - WS * (1 - us))
  RS = S + WS * us
  RE = pmin(RS, at + WE * ue)
  data.frame(LE,LS,RE,RS)  
}

trueparam <- list(shape=1.5, scale=10)
sigma <- 1.0
tau <- 100L
set.seed(1111)
curve(Lambda(x,sigma),0,100)
at = weibull_process(tau, shape = sigma, maxit = 10000)
lines(at[-1], 1:length(at[-1]), type = "s", col="grey")

y <- rweibull(length(at), trueparam$shape, scale = trueparam$scale)

WE=1
WS=1


simDat0 <- simDIC(y, at = at, WS = 1, WE = 1)
simDat <- simDat0[simDat0$RS < tau,]

system.time({
  est <- jointfit_dic(LE = simDat$LE, RE = simDat$RE, LS = simDat$LS,　RS =  simDat$RS, tau = tau,
                      maxit=100L, tol=1e-8, method = "VB")
})

which(is.nan(est$alpha))
which(is.nan(est$beta))
est$logprob
plot(est$breaks, cumsum(est$incubation), type="s")



system.time({
  est2 <- jointfit_dic(LE = simDat$LE, RE = simDat$RE, LS = simDat$LS,　RS =  simDat$RS, tau = tau,
                      maxit=100L, tol=1e-8, method = "EM")
})

plot(est2$breaks, cumsum(est2$incubation), type="s")
curve(pweibull(x, shape = trueparam$shape, scale = trueparam$scale), add=TRUE, lty=2)

system.time({
  res <- joint_dic_gibbs(LE = simDat$LE, RE = simDat$RE, LS = simDat$LS,　RS =  simDat$RS, tau = tau,
                       iter=1000L)
})

plot(res$logprob, type = "l")

dim(res$incubation)

burnin <- 1:100
Shat = rowMeans(1-apply(res$incubation[,-burnin],2,cumsum))
plot(res$breaks, Shat, type="s", col="royalblue", xlim = c(0,100))
curve(pweibull(x, shape = trueparam$shape, scale = trueparam$scale, lower.tail = FALSE), add=TRUE, lty=2)


Lamhat = rowMeans(res$infection[,51:100])
plot(res$breaks, Lamhat, type="s", col="royalblue", xlim = c(0,100))
abline(0,1)

plot(est$breaks, 1-cumsum(est$incubation), type="s", col="royalblue", xlim = c(0,30))
lines(est2$breaks, 1-cumsum(est2$incubation), type="s", col="orangered")
lines(res$breaks, Shat, type="s", col="forestgreen")
curve(pweibull(x, shape = trueparam$shape, scale = trueparam$scale, lower.tail = FALSE), add=TRUE, lty=2)

####
df_post <- IIdf(est)
ggplot(df_post, aes(x=time))+
  geom_step(aes(y=cuminf, colour="not tuncated"))+
  geom_step(aes(y=cuminf_unobs, colour="tuncated"))+
  theme_bw(16)+
  labs(fill="",colour="",y="cumulative infections")

set.seed(1);rand = rand_cdf(est$beta, np = 2000L)
dfrand_i <- melt(rand, varnames = c("index","group"), value.name = "intensity") %>% 
  mutate(x=est$breaks[index])
df_i <- data.frame(x = est$breaks, intensity=cumsum(est$beta))

p1 <- ggplot(df_i, aes(x=x, y = intensity))+
  geom_step()+
  geom_step(data=dfrand_i, aes(group=group), alpha=0.01)+
  stat_function(fun = function(x){x^sigma}, colour="orangered", linetype=2)+
  #geom_abline(slope = 1, intercept = 0, colour="orangered", linetype=2)+
  theme_bw(15)+labs(title = "infection")

print(p1)

dfrand_in <- melt(normalize_cdf(rand), varnames = c("index","group"), value.name = "intensity") %>% 
  mutate(x=est$breaks[index])
df_in <- data.frame(x = est$breaks, intensity=cumsum(est$beta/sum(est$beta)))

p2 <- ggplot(df_in, aes(x=x, y = intensity))+
  geom_step()+
  geom_step(data=dfrand_in, aes(group=group), alpha=0.01)+
  stat_function(fun = function(x){(x^sigma)/Lambda(tau,sigma)}, colour="orangered", linetype=2)+
  #geom_abline(slope = 1/n, intercept = 0, colour="orangered", linetype=2)+
  theme_bw(15)+labs(title = "infection (relative)", y="intensity (relative)")

print(p2)


set.seed(2); rand2 = rand_cdf(est$alpha, np = 2000L)
dim(rand2)
dfrand2_i
dfrand2_i <- melt(rand2, varnames = c("index","group"), value.name = "intensity") %>% 
  mutate(x=est$breaks[index])
df2_i <- data.frame(x = est$breaks, intensity=cumsum(est$alpha))


p3 <- ggplot(df2_i, aes(x=x, y = intensity))+
  geom_step()+
  geom_step(data=dfrand2_i, aes(group=group), alpha=0.01)+
  stat_function(fun = function(x)Lambda(tau, sigma)*pgamma(x, shape = trueparam$shape, scale = trueparam$scale), colour="royalblue", linetype=2)+
  #stat_function(fun = function(x)Lambda(tau, sigma)*pweibull(x, shape = trueparam$shape, scale = trueparam$scale), colour="royalblue", linetype=2)+
  theme_bw(15)+labs(title = "incubasion")
print(p3)

dfrand2_in <- melt(normalize_cdf(rand2), varnames = c("index","group"), value.name = "intensity") %>% 
  mutate(x=est$breaks[index])
df2_in <- data.frame(x = est$breaks, intensity=cumsum(est$incubation))

p4 <- ggplot(df2_in, aes(x=x, y = intensity))+
  geom_step()+
  geom_step(data=dfrand2_in, aes(group=group), alpha=0.01)+
  stat_function(fun = function(x)pgamma(x, shape = trueparam$shape, scale = trueparam$scale), colour="royalblue", linetype=2)+
  #stat_function(fun = function(x)pweibull(x,shape = trueparam$shape, scale = trueparam$scale), colour="royalblue", linetype=2)+
  theme_bw(15)+labs(title = "incubasion (relative)", y="intensity (relative)")

print(p4)


# pc <- grid.arrange(p1, p2, p3, p4)
# ggsave("rand_int.pdf", plot = pc)


tp = c(0.5,0.8,0.9)
brks <- qweibull(tp, trueparam$shape, trueparam$scale)
tvalue = c(0.5,0.8,0.9)*50
#50*pweibull(brks,2,5)
df <- data.frame(x = est$breaks, 
                 estimates = apply(rand2, 1, mean))

dfci <- simulated_pvalue(rand2, est, tp, brks, tvalue)
ggplot(df, aes(x=x, y=estimates))+
  stat_function(fun = function(x){n*pweibull(x,trueparam$shape,trueparam$scale)}, colour="orangered", linetype=2)+
  geom_pointrange(data = dfci, aes(ymin=lower, ymax=upper, colour=pvalue))+
  geom_step()+
  scale_color_continuous(type = "viridis")+
  theme_bw()

ggsave("pvalue.pdf")
save.image()

