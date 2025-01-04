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
#save.image()
plot(est$logprob, type="l")

k01 = rgb(0,0,0,0.01)

plot(est$breaks,cumsum(est$beta), type="s", col="royalblue", ylim=c(0,100),xlim=c(0,100))
abline(0,1,col="orangered", lty=2)
rand = rand_cdf(est$beta, np = 1000L)
matlines(est$breaks, rand, col=k01, lty=1)

#png("ccdf.png")
#par(mfrow=c(2,1))
curve(50*pweibull(x,2,5), type="l", col="orangered", xlim = range(est$breaks),ylim=c(0,80))
rand = rand_cdf(est$alpha, np = 1000)
matlines(est$breaks, rand, col=k01, lty=1)
curve(pweibull(x,2,5), type="l", col="orangered", xlim = range(est$breaks))
matlines(est$breaks, normalize_cdf(rand), type="s", lty = 1,col=k01)
#dev.off()

simulated_pvalue <-function(rand, est, tp, brks, tvalue){
  CI <- matrix(0, nrow = length(brks), ncol = 3)
  colnames(CI) <- c("lower", "upper", "estimates")
  
  xv = numeric(length(brks))
  pv <- numeric(length = length(brks))
  for(i in 1:length(brks)){
    ind = which.max(est$breaks > brks[i])-1
    ps <- seq(0, 1, length.out=ncol(rand))
    rs <- sort(rand[ind,])
    lower <- ps[which.min(rs < tvalue[i])]
    upper <- rev(ps)[which.min(tvalue[i]<rev(rs))]
    pv[i] <- min(lower, upper)
    
    xv[i] <- est$breaks[ind]
    ci <- quantile(rand[ind,], probs = c(pv[i], 1-pv[i], 0.5))
    CI[i,] <- ci
  }
  data.frame(x=brks, pvalue=pv, CI)
}

df <- data.frame(x = est$breaks, 
                 estimates = apply(rand, 1, mean))

tp = c(0.5,0.8,0.9)
brks <- qweibull(tp, 2,5)
tvalue = c(0.5,0.8,0.9)*50
#50*pweibull(brks,2,5)

dfci <- simulated_pvalue(rand, est, tp, brks, tvalue)
ggplot(df, aes(x=x, y=estimates))+
  stat_function(fun = function(x){50*pweibull(x,2,5)}, colour="orangered", linetype=2)+
  geom_pointrange(data = dfci, aes(ymin=lower, ymax=upper, colour=pvalue))+
  geom_step()+
  scale_color_continuous(limits=c(0, 1))+
  theme_bw()

#ggsave("pvalue.pdf")
  