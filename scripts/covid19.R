library(readr)
library(dplyr)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(JuliaCall)
library(ggplot2)
library(JointII)
path <- scan("path.txt", what = character())
julia <- julia_setup(JULIA_HOME = path[1])
julia_source(path[2])
julia_source(path[3])

#####
#pre treat
dat_raw <- read_csv(path[4], skip = 1)
head(dat_raw)
dat_out <- dat_raw %>%
  dplyr::filter(!is.na(exposure_end)) %>% 
  mutate(exposure_start=as.Date(exposure_start,format = "%m/%d/%y"),
         exposure_end=as.Date(exposure_end,format = "%m/%d/%y"),
         symptom_onset=as.Date(symptom_onset,format = "%m/%d/%y"),
         hosp_visit_date=as.Date(hosp_visit_date,format = "%m/%d/%y"),
         reporting_date=as.Date(`reporting date`,format = "%m/%d/%y")) %>%
  dplyr::select(id,exposure_start,exposure_end,symptom_onset,hosp_visit_date,reporting_date) %>% 
  dplyr::filter(ifelse(is.na(exposure_start),TRUE,exposure_start<=exposure_end)) %>% 
  mutate(LE = exposure_start,
         RE = exposure_end,
         LS = if_else(!is.na(symptom_onset), symptom_onset, LE),
         RS = if_else(!is.na(symptom_onset), symptom_onset, if_else(!is.na(hosp_visit_date), hosp_visit_date, reporting_date)))


dat4 <- dplyr::select(dat_out,LE,RE,LS,RS) %>% 
  dplyr::filter(!is.na(LE)) %>% 
  mutate(recmin = min(c(LE,RE,LS,RS))) %>% 
  mutate(LE=as.integer(LE-recmin), 
         RE=as.integer(RE-recmin),
         LS=as.integer(LS-recmin), 
         RS=as.integer(RS-recmin)) %>% 
  dplyr::filter(!((RE>=LS & LS==RS) | (LE >= LS & LE==RE))) %>% 
  mutate(ctype=if_else(LS==RS & LE==RE, 0L, 3L)) %>% 
  mutate(ctype=if_else(LS==RS & LE<RE, 1L, ctype)) %>% 
  mutate(ctype=if_else(LS<RS & LE==RE, 2L, ctype)) 

tau <- max(dat4$RS)+1

julia_command("using Random; Random.seed!(1234)")
res2 <- julia_call("jointest.gibbssampler", dat4$LE, dat4$RE, dat4$LS, dat4$RS, tau, 2000L)
plot(res2$logprob[-1], type = "l")

system.time({
  res <- joint_dic_gibbs(LE = dat4$LE, RE = dat4$RE, LS = dat4$LS,　RS =  dat4$RS, tau = tau,
                         iter=2000L)
})

burnin <- 1:1000
Shat = rowMeans(1-apply(res$incubation[,-burnin],2,cumsum))

plot(colMeans(res2$ccdf[-burnin,]), type = "s")
lines(Shat, type = "s", col="blue")

df_eccdf <- data.frame(
  value = res2$value,
  ccdf = colMeans(res2$ccdf[1001:2000, ])
)

res2 <- julia_call("empest.gibbssampler", dat4$LE, dat4$RE, dat4$LS, dat4$RS, 2000L)

ggplot(df_eccdf, aes(x=value, y=ccdf))+
  geom_step()+
  theme_classic()
