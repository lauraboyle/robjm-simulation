 
# simulate normal normal mixed model code 
# similar to size of renal 

setwd("~/Documents/GitHub/robjm-simulation")
source("Simulate_Data_nor_nor.R")
source("Supplementary_functions.R")


gamma<-0
nsubj<-600
av_n_i <- 12

simulated_data <- simulate_data(nsubj = nsubj, av_n_i = av_n_i, gamma=gamma, returns = c("repeat_data", "base_data"))

repeat_data<-simulated_data$repeat_data
dim(repeat_data)

fit_nor_nor <- fit_ld(fixed = Y ~ time , 
                      random = ~ time , 
                      data = repeat_data,
                      id = "id",
                      model = "nor_nor",
                      chains = 2, 
                      cores=2,
                      iter = 2000,
                      warmup = 1000,
                      control = list(adapt_delta = 0.99))
print(fit_nor_nor, pars = c("alpha", "Sigma", "sigmasq"))