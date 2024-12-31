#weibull_process(100, shape = 1.5, maxit = 1000)
y <- rweibull(100, 2, 5)
simDat <- simDIC(y)
head(simDat)
# breaks <- setbreaks(simDat$LE, simDat$RE, simDat$LS, simDat$RS)
# RANK <- rankbreaks(simDat$LE, simDat$RE, simDat$LS, simDat$RS, breaks)
est <- jointecdf_dic_em(simDat$LE, simDat$RE, simDat$LS, simDat$RS,
                        maxit=20L, tol=1e-5)
plot(est$lp)
# RANK$LE[1]:RANK$RE[1]
# RANK$LS[1]:RANK$RS[1]
plot(est$breaks, p2ccdf(est$alpha/sum(est$alpha)), type = "s")
curve(pweibull(x,2,5, lower.tail = FALSE), add=TRUE,lty=2)
plot(est$breaks, cumsum(est$beta), type = "s")
abline(0,1, lty=2)
image(est$breaks, est$breaks, est$event, col = grey.colors(12))
