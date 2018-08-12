
 # loading the packages to be used

 library(rstudioapi) 
 library(ggplot2)
 library(JM)
 library(tidyverse)
 library(mvtnorm)
 library(Matrix)  
 library(survminer)

 # set the wd to where the codes are located
 setwd(dirname(getActiveDocumentContext()$path))
 
 # source the supp. functions
 source("Supplementary_functions.R")
 source("Simulate_Data_nor_nor.R")
 source("Simulate_Data_tt_mod3.R")
 source("Simulate_Data_tt_tv.R")
 
 # number of subjects for modelling and external validation
 nsubj  <- 100
 phi    <- 4
 delta0 <- 6
 beta <- c(-2, 1, 3)
 gamma <- 0.3
 
 # simulate data to be used for parameter estimation
# simulated_data <- simulate_data(nsubj = nsubj, av_n_i = av_n_i, returns = c("repeat_data", "base_data"))
# simulated_data_tt_mod3 <- simulate_data_tt_mod3(nsubj = nsubj, av_n_i = av_n_i, phi=phi, delta=delta, returns = c("repeat_data", "base_data"))
 simulated_data <- simulate_data_tt_tv(nsubj = nsubj, 
                                       av_n_i = av_n_i, 
                                       phi = phi, 
                                       gamma = gamma,
                                       delta0 = delta0, 
                                       beta = beta,
                                       returns = c("repeat_data", "base_data"))
 
 # mixed-model
 lme_fit <- lme(fixed = Y ~ time, random = ~ time|id, data = simulated_data$repeat_data)
 summary(lme_fit)
 
 cox_fit <- coxph(Surv(stime, event) ~ c, data = simulated_data$base_data, x = T)
 summary(cox_fit)
 
 joint_fit <- jointModel(lme_fit, cox_fit, timeVar = "time",method = "weibull-PH-aGH",
                         verbose = T, iter.EM = 1000)
 summary(joint_fit)
 
 g <- ggplot(simulated_data$repeat_data, aes(time, Y, group = id))
 g + geom_point()
 g + geom_line()
  
 summary(unlist(with(simulated_data$repeat_data, tapply(time, id, diff))))
  
 summary(as.numeric(table(simulated_data$repeat_data$id)))
 
 km_fit <- survfit(Surv(stime, event) ~ 1, data = simulated_data$base_data)
 ggsurvplot(km_fit, data = simulated_data$base_data, risk.table = TRUE)
 
 
 
 # fit using robjm 
 library(robjm)
 fit_nor_nor <- fit_jm(fixed_long = Y ~ time, 
                       random_long = ~ time, 
                       fixed_surv = cbind(stime, event)~c, 
                       data_long = simulated_data$repeat_data,
                       data_surv = simulated_data$base_data,
                       id_long = "id",
                       id_surv = "id",
                       model = "nor_nor",
                       chains = 2,
                       cores=2,
                       iter = 2000,
                       warmup = 1000,
                       control = list(adapt_delta = 0.99)
 )
 print(fit_nor_nor, pars = c("alpha", "Sigma", "sigmasq", "zeta", "omega", "eta"))

 fit_t_t_tv <- fit_jm(fixed_long = Y ~ time, 
                      random_long = ~ time,  
                      fixed_surv = cbind(stime, event)~c, 
                      data_long = simulated_data$repeat_data,
                      data_surv = simulated_data$base_data,
                      id_long = "id",
                      id_surv = "id",
                      model = "t_t_tv",
                      spline_tv = list("time", 2), 
                      chains = 2,
                      iter = 2000,
                      warmup = 1000,
                      control = list(adapt_delta = 0.99)
 )
 print(fit_t_t_tv, pars = c("alpha", "Sigma", "phi", "sigmasq", "delta0", "beta", "zeta", "omega", "eta")) 
 