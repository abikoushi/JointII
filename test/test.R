library(JointII)
library(Matrix)
library(bench)
library(ggplot2)
library(pammtools)
library(reshape2)
library(gridExtra)
#library(RColorBrewer)

logsumexp <- function(x){
  mx <- max(x)
  mx + log(sum(exp(x-mx)))
}

logsumexp(1:10)
log(sum(exp(1:10)))

Lambda <- function(x, sigma){
  x^sigma
}

trueparam <- list(shape=1.5, scale=3)
sigma <- 1.1
tau <- 100L
set.seed(7689)
at = weibull_process(tau, shape = sigma, maxit = 10000)
y <- rgamma(length(at), trueparam$shape, scale = trueparam$scale)
simDat0 <- simDIC(y, at = at)
simDat <- simDat0[simDat0$RS<tau, ]

system.time({
  est <- jointfit_dic(simDat$LE, simDat$RE, simDat$LS, simDat$RS, tau = tau,
                          maxit=2L, tol=0, method = "VB")
})

system.time({
  est2 <- jointfit_dic(simDat$LE, simDat$RE, simDat$LS, simDat$RS, tau = tau,
                      maxit=3L, tol=0, method = "VB")
})

plot(est$logprob)
image(t(est$truncated_event))

range(est$alpha[!is.nan(est2$alpha)])
range(est$alpha[is.nan(est2$alpha)])

range(est$beta[!is.nan(est2$beta)])
range(est$beta[is.nan(est2$beta)])

plot(est$breaks, cumsum(est$incubation), type="s")
curve(pgamma(x, shape = trueparam$shape, scale = trueparam$scale), col="royalblue", add=TRUE)

plot(est2$breaks, cumsum(est2$incubation), type="s")
curve(pgamma(x, shape = trueparam$shape, scale = trueparam$scale), col="royalblue", add=TRUE)


plot(est$logprob, type="l")
dim(simDat0)
dim(simDat)

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

