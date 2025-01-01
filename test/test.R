library(JointII)
library(Matrix)
library(bench)
#library(RColorBrewer)
library(ggplot2)
library(pammtools)
#weibull_process(100, shape = 1.5, maxit = 1000)
y <- rweibull(50, 2, 5)
simDat <- simDIC(y)
head(simDat)
# breaks <- setbreaks(simDat$LE, simDat$RE, simDat$LS, simDat$RS)
# RANK <- rankbreaks(simDat$LE, simDat$RE, simDat$LS, simDat$RS, breaks)
bmrec <- bench::mark({
  est <- jointecdf_dic_em(simDat$LE, simDat$RE, simDat$LS, simDat$RS,
                          maxit=200L, tol=0.01)
}, min_iterations = 1, max_iterations = 1)

plot(est$logprob, type="l")

curve(pweibull(x,2,5), type="l", col="orangered", xlim = range(est$breaks))
lines(est$breaks,cumsum(est$alpha)/sum(est$alpha), type="s", col="royalblue")
out = numeric(length(est$alpha))
out[wch] = rgamma(length(wch), est$alpha[wch], 1)
out = out/sum(out)
lines(est$breaks,cumsum(out), type="s", col=rgb(0,0,0,0.1))

plot(est$logprob, type = "l")
# bmrec2 <- bench::mark({
#   est2 <- jointecdf_dic_vb(simDat$LE, simDat$RE, simDat$LS, simDat$RS,
#                           maxit=1000L, tol=1e-5)    
# }, min_iterations = 1, max_iterations = 1)

wch <- which(est$incubation > 0)
m <- length(est$incubation[wch])
I <- matrix(1,m,m)
diag(I) <- 2/est$incubation[wch]-1
I <- sum(est$event)*I
vhat <- numeric(length = length(est$incubation))
vhat[wch] <- diag(I)

se <- sqrt(p2ccdf(vhat))

length(se)
dim(resdf)
resdf <- IIdf(est)
tail(resdf)
ggplot(resdf, aes(x=time, y=ccdf))+
  geom_step()+
  #geom_stepribbon(aes(ymin=ccdf-se, ymax=ccdf+se), alpha=0.1)+
  stat_function(fun = pweibull, args = list(shape=2, scale=5, lower.tail = FALSE), linetype=2, colour="royalblue")+
  theme_bw()

bmrec2$total_time

image(est$event + est$truncated_event)

plot(est$logprob, type = "l")

resdf <- IIdf(est)
ggplot(resdf, aes(x=time, y=ccdf))+
  geom_step()+
  stat_function(fun = pweibull, args = list(shape=2, scale=5, lower.tail = FALSE), linetype=2, colour="royalblue")+
  theme_bw()

ggplot(resdf, aes(x=time, y=intensity))+
  geom_step()+
  stat_function(fun = function(x)x, linetype=2, colour="orangered")+
  theme_bw()

ggplot(resdf, aes(x=time, ))+
  geom_step(aes(y=cuminf+cuminf_unobs))+
  geom_step(aes(y=cuminf_unobs), colour="grey")+
  theme_bw()

# plot(resdf$time, resdf$ccdf, type="s")
# curve(pweibull(x, 2, 5, lower.tail = FALSE), add=TRUE, col="royalblue", lty=2)

# plot(resdf$time, resdf$intensity, type="s")
# abline(0, 1, col="royalblue", lty=2)

# plot(resdf$time, resdf$cuminf, type="s")
# lines(resdf$time, resdf$cuminf_unobs, type="s", col="orangered", lty=2)
