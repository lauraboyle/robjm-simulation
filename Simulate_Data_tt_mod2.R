simulate_data_tt_mod2 <- function(nsubj , av_n_i = 5, gamma , phi , delta,  returns = c("repeat_data", "base_data")){
  
  ### nsubj: number of subjects
  ### delta: the parameter that controls informative visit process
  ### return: names of the objects to be returned
  
  # discretise time at fine intervals
  t_min <- 0
  t_max <- 5
  a     <- 0.01
  t     <- seq(t_min, t_max, a)
  m     <- length(t)
  
  # total number of observations before censoring and selection
  ntotal <- m * nsubj
  
  # longitudinal fixed effects regression parameters
  alpha <- matrix(c(1, 0.5), ncol = 1)
  
  # var-cov matrix of random intercept and slope
  Sigma <- matrix(c(0.6, 0.25, 0.25, 0.3), ncol = 2) 
  
  # measurement error standard deviation
  tau <- sqrt(0.25)
  
  # longitudinal fixed-effects covariate matrix
  x <- cbind(1, rep(t, nsubj)) 
  
  # random effects covariate matrix
  d <- Matrix::bdiag(lapply(1:nsubj, function(i) cbind(1, t)))
  
  # random effects
  U_mat <- mvtnorm::rmvnorm(nsubj, mean = rep(0, ncol(Sigma)), sigma = Sigma)
  U     <- as.numeric(t(U_mat))
  
  # measurement error
  Z <- rnorm(ntotal, 0, tau)
  
  # create gamma distributed Vi                             
  V <- rgamma(nsubj, shape = phi/2, rate= phi/2)              
  V_ext <- rep(V,each=501)                                    
  
  # create gamma distributed Wi                              
  W <- rgamma(nsubj, shape = delta/2, rate= delta/2)              
  W_ext <- rep(W,each=501)
    
  # longitudinal measurements w/o error
  #  Y_star <- as.matrix(x %*% alpha + d %*% U)               
  Y_star <- as.matrix(x %*% alpha + ((d %*% U)/sqrt(V_ext)))
  
  # longitudinal measurements with error
  #  Y <- Y_star + Z                                          
  Y <- Y_star + (Z/sqrt(W_ext))                                   
  
  # data-set
  data <- data.frame(
    id   = rep(1:nsubj, each = m),
    time = x[, 2],
    Y = Y,
    Y_star = Y_star,
    U0 = rep(U_mat[, 1], each = m),
    U1 = rep(U_mat[, 2], each = m)
  )
  
  # survival process parameters
  theta <- 0.04
  phi   <- 1.2
  omega <- matrix(0, ncol = 1, nrow = 1)
  #gamma <- 0.3
  
  # survival model covariates
  c <- matrix(rep(0, ntotal))
  
  # hazard and survival probabilities
  data$hazard <- as.numeric(with(data, theta * phi * time^(phi - 1) * exp(c %*% omega + gamma * Y_star)))
  a_vec <- c(1, rep(a, (m - 1)))
  data$surv_prob <- unlist(with(data, tapply(hazard, id, function(x) exp(- unlist(lapply(1:m, function(i) sum(a_vec[1:i] * x[1:i])))))))
  
  # select observations randomly by being sure that everyone has data at baseline, 0: non-selection, 1:selection
  data$sel <- unlist(lapply(1:nsubj, function(i) c(1, rbinom((m - 1), 1, (av_n_i - 1)/(m - 1)))))
  
  # uniform random variables for survival 
  data$unif_rand_survival  <- rep(runif(nsubj), each = m)
  
  # censor at the event
  data_censored_event <- data[data$surv_prob > data$unif_rand_survival, ]
  
  # create stime and event indicator
  stime_event <- 
    do.call(rbind, with(data_censored_event, tapply(time, id, function(x) stime_event_fun(x, t_max = t_max, a = a))))
  data_censored_event$stime <- stime_event[, 1]
  data_censored_event$event <- stime_event[, 2]
  
  # further select visits informatively
  data_censored_event_sel <- dplyr::filter(data_censored_event, sel == 1)
  
  # final data-set: repeat + base 
  repeat_data <- data_censored_event_sel
  base_data   <- repeat_data[repeat_data$time == 0, ]
  
  # return some data
  mget(returns)
  
}
