library(ggplot2)
library(Matrix)
library(JointII)
library(JuliaCall)
juliapath <- scan("path.txt", what = character())
julia <- julia_setup(JULIA_HOME = juliapath[1])
julia_source(juliapath[2])
julia_source(juliapath[3])

#累積強度関数
Intensity <- function(t,a,b){
  (exp(a*t+b)-exp(b))/a
}

#強度関数
intensity <- function(t,a,b){
  exp(a*t+b)
}

NHPP_exp = function(Tmax, a, b, maxit){
  mlogz = -log(runif(maxit))
  t <- rep(Inf, maxit)
  t[1] <- log(a*mlogz[1]-b + exp(b))/a
  for(i in 2:maxit){
    s <- a*mlogz[i] + exp(a*t[i-1]+b) 
    if(s<0){
      break
    }
    ti =(log(s) - b)/a
    if(ti>Tmax){
      break
    }
    t[i] <- ti
  }
  return(t[1:(i-1)])
}


ti <- NHPP_exp(100, 0.01, 0, maxit=1000)
y <- rgamma(length(ti), shape = 2, scale = 5)
S <- ti+y
obs <- S<100
df <- data.frame(E = ti[obs], S = S[obs])
system.time({
  est_vb <- jointfit_dic(LE = floor(df$E), RE = ceiling(df$E),
                         LS = floor(df$S),　RS = ceiling(df$S),
                         tau = 100,
                         maxit=1000L, tol=1e-8, method = "VB")
})

plot(est_vb$breaks, cumsum(est_vb$beta))
curve(Intensity(x,0.01,0), add=TRUE)

est_vb2 <- julia_call("jointtool.vb",
                      floor(df$E), ceiling(df$E),
                      floor(df$S), ceiling(df$S), 100.0)

# est_vb2 <- julia_call("jointtool.vb", 
#                       df$E, df$E,
#                       df$S, df$S, 100.0)

plot(est_vb2$value, cumsum(est_vb2$beta))
curve(Intensity(x,0.01,0), add=TRUE)

plot(est_vb2$value, cumsum(est_vb2$alpha/sum(est_vb2$alpha)))
curve(pgamma(x, shape = 2, scale = 5), add=TRUE)

points(est_vb$breaks, cumsum(est_vb$incubation))
curve(pgamma(x, shape = 2, scale = 5), add=TRUE)

plot(est_vb$breaks, cumsum(est_vb$beta))
curve(Intensity(x,0.01,0), add=TRUE)

df_est <- data.frame(breaks = est_vb$breaks,
                 B=cumsum(rowSums(est_vb$truncated_event)),
                 A=cumsum(rowSums(est_vb$event)))

ggplot(df_est, aes(x=breaks))+
  geom_step(aes(y=A+B))+
  geom_step(aes(y=B, colour="B"))+
  theme_bw()

