
# inverse of logit
expit <- function(x){
  exp(x)/(1+exp(x))
}

# function to create surv time and event indicator
stime_event_fun <- function(x, t_max, a){
  if(max(x) < t_max){
    matrix(c(max(x) + a, 1), ncol = 2)[rep(1, length(x)), ]
  }else{
    matrix(c(max(x), 0), ncol = 2)[rep(1, length(x)), ]
  }
}
