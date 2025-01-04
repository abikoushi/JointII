library(JointII)
library(Matrix)
library(bench)
library(ggplot2)
library(pammtools)
library(reshape2)
library(gridExtra)
#library(RColorBrewer)

#weibull_process(100, shape = 1.5, maxit = 1000)

trueparam <- list(shape=2, scale=5)
n <- 50
set.seed(333)
y <- rweibull(n, trueparam$shape, trueparam$scale)
simDat <- simDIC(y)

bmrec <- bench::mark({
  est <- jointfit_dic(simDat$LE, simDat$RE, simDat$LS, simDat$RS,
                          maxit=1000L, tol=1e-5, method = "VB")
}, min_iterations = 1, max_iterations = 1)

plot(est$logprob, type="l")


set.seed(1);rand = rand_cdf(est$beta, np = 1000L)
dfrand_i <- melt(rand, varnames = c("index","group"), value.name = "intensity") %>% 
  mutate(x=est$breaks[index])
df_i <- data.frame(x = est$breaks, intensity=cumsum(est$beta))

p1 <- ggplot(df_i, aes(x=x, y = intensity))+
  geom_step()+
  geom_step(data=dfrand_i, aes(group=group), alpha=0.01)+
  geom_abline(slope = 1, intercept = 0, colour="orangered", linetype=2)+
  theme_bw(15)+labs(title = "infection")

print(p1)

dfrand_in <- melt(normalize_cdf(rand), varnames = c("index","group"), value.name = "intensity") %>% 
  mutate(x=est$breaks[index])
df_in <- data.frame(x = est$breaks, intensity=cumsum(est$beta/sum(est$beta)))

p2 <- ggplot(df_in, aes(x=x, y = intensity))+
  geom_step()+
  geom_step(data=dfrand_in, aes(group=group), alpha=0.01)+
  geom_abline(slope = 1/n, intercept = 0, colour="orangered", linetype=2)+
  theme_bw(15)+labs(title = "infection (relative)", y="intensity (relative)")

print(p2)

set.seed(2); rand2 = rand_cdf(est$alpha, np = 1000L)
dfrand2_i <- melt(rand2, varnames = c("index","group"), value.name = "intensity") %>% 
  mutate(x=est$breaks[index])
df2_i <- data.frame(x = est$breaks, intensity=cumsum(est$alpha))

p3 <- ggplot(df2_i, aes(x=x, y = intensity))+
  geom_step()+
  geom_step(data=dfrand2_i, aes(group=group), alpha=0.01)+
  stat_function(fun = function(x)n*pweibull(x,shape = trueparam$shape, scale = trueparam$scale), colour="royalblue", linetype=2)+
  theme_bw(15)+labs(title = "incubasion")

print(p3)

dfrand2_in <- melt(normalize_cdf(rand2), varnames = c("index","group"), value.name = "intensity") %>% 
  mutate(x=est$breaks[index])
df2_in <- data.frame(x = est$breaks, intensity=cumsum(est$incubation))

p4 <- ggplot(df2_in, aes(x=x, y = intensity))+
  geom_step()+
  geom_step(data=dfrand2_in, aes(group=group), alpha=0.01)+
  stat_function(fun = function(x)pweibull(x,shape = trueparam$shape, scale = trueparam$scale), colour="royalblue", linetype=2)+
  theme_bw(15)+labs(title = "incubasion (relative)", y="intensity (relative)")

print(p4)

pc <- grid.arrange(p1, p2, p3, p4)
ggsave("rand_int.pdf", plot = pc)


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

