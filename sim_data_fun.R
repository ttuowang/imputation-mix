#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Title:  Imputation of a mix of variable types                            #
# Description: This program contains data generation procedure             #    
# Author: Tuo Wang (email: twang437@wisc.edu, website: tuowang.rbind.io)   #
# Date:   July-Aug 2021                                                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

sim_event <- function(
  A, X, b1, b2, 
  lambda=1, k1 , k2, alpha, beta1, beta2, 
  recurrent=FALSE, tmax =12, start=0, end=1000){
  
  # This is the function for data generation for one subject
  
  if(recurrent == FALSE){
    Delta <- lambda * exp(
      alpha[3]*beta1[1] + alpha[4]*beta2[1] + alpha[1]*A +
        (alpha[2] + alpha[3]*beta1[2] + alpha[4]*beta2[2]) * X
    )
    u <- runif(1, 0, 1)
    f1 <- function(s){
      exp(
        alpha[3]*(beta1[3]*A+b1)*(1 - exp(-k1*s)) + 
          alpha[4]*(beta2[3]*A+b2)*(1 - exp(-k2*s)) )
    }
    f2 <- function(t){
      int <- integrate(f1, lower = 0, upper = t)
      (-log(u))/(Delta) - int$value
    }
    tryCatch(
      {TD <- uniroot(f2, c(start, end))$root},
      error = function(err){
        print(paste(err, "Terminal - Try increase the end point"))}
    )
    return(TD)
  }
  if(recurrent == TRUE){
    TH <- 0
    while (tail(TH, 1) < tmax) {
      Delta <- lambda * exp(
        alpha[3]*beta1[1] + alpha[4]*beta2[1] + alpha[1]*A +
          (alpha[2] + alpha[3]*beta1[2] + alpha[4]*beta2[2]) * X
      )
      u <- runif(1, 0, 1)
      f1 <- function(s){
        exp(
          alpha[3]*(beta1[3]*A+b1)*(1 - exp(-k1*s)) + 
            alpha[4]*(beta2[3]*A+b2)*(1 - exp(-k2*s)) )
      }
      f2 <- function(t){
        int <- integrate(f1, lower = 0, upper = t)
        (-log(u))/(Delta) - int$value
      }
      tryCatch(
        {W <- uniroot(f2, c(start, end))$root},
        error = function(err){
          print(paste(err, "Recurrent - Try increase the end point"))}
      )
      if(W > 1e-3){
        if( tail(TH,1) == 0){
          if(tail(TH,1)+W > tmax){TH = NULL; break}
          TH <- W
        }else{
          if(tail(TH,1)+W > tmax){break}
          TH <- c(TH, tail(TH,1)+W)
        }
      }
    }
  }
  return(TH)
}


sim_data <- function(
  n, tmax, t, 
  k1, k2, beta1, beta2, 
  alpha, gamma, lambdaD0=1, lambdaH0=1 , 
  recurrent=TRUE, censor="independent",
  start=0, end=1000, seed=NULL){
  
  # This is the function for data generation for multiple subjects
  # n: sample size
  # tmax: maximum follow-up time
  # t: time points to measure longitudinal data
  # k1, k2: time effect for cts and binary var respectively
  # beta1, beta2: fixed effect for cts and binary var
  # alpha: coefficients for time-to-event
  # gamma: coefficients for recurrent events
  # lambdaD0: baseline hazard for time-to-event
  # lambdaH0: baseline hazard for recurren events
  
  expit <- function(x){
    exp(x)/(1+exp(x))
  }
  
  set.seed(seed)
  
  A <- rbinom(n, 1, 0.5)
  X <- round(rnorm(n, 0, 1),4)
  
  mu.b <- c(0, 0)
  sigma2.b1 <- 0.1
  sigma2.b2 <- 0.1
  rho <- 0.5
  cov.b <- rho * sqrt(sigma2.b1 * sigma2.b2)
  cov.mat.b <- matrix(c(sigma2.b1, cov.b, cov.b, sigma2.b2), nrow=2)
  
  b.mat <- MASS::mvrnorm(n, mu.b, cov.mat.b)
  b1 <- b.mat[, 1]
  b2 <- b.mat[, 2]
  
  nrm <- length(t)   
  
  ft1 <- outer((beta1[3]*A + b1), (1- exp(-k1*t)) , FUN = "*")
  mt1 <- beta1[1] + beta1[2]*X + ft1
  epsilon1 <- matrix(rnorm(n*nrm, mean = 0, sd = 0.4), nrow = n, ncol = nrm)
  Y1 <- mt1 + epsilon1
  Y1 <- round(Y1,4)
  
  ft2 <- outer((beta2[3]*A + b2), (1- exp(-k2*t)) , FUN = "*")
  mt2 <- beta2[1] + beta2[2]*X + ft2
  epsilon2 <- matrix(rnorm(n*nrm, mean = 0, sd = 0.4), nrow = n, ncol = nrm)
  Y2 <- matrix( 
    1*(expit(mt2 + epsilon2) > 0.5), 
    nrow = n, ncol = nrm)

  cov_data <- data.frame(
    id = c(1:n),
    A = A,
    X = X,
    Y1,
    Y2
  )
  names(cov_data) <- c(
    "id", "A", "X", paste0("Y1_t",t), paste0("Y2_t",t)
  )
  
  event_data <- NULL
  event_miss <- NULL
  for (i in 1:n) {
    
    Ai <- A[i]
    Xi <- X[i]
    b1i <- b1[i]
    b2i <- b2[i]

    # Terminal events
    TD <- sim_event(
      A=Ai, X=Xi, b1=b1i, b2=b2i, lambda=lambdaD0,
      k1=k1, k2=k2, alpha=alpha, beta1=beta1, beta2=beta2,
      recurrent=FALSE, start=start, end=end
    )
    
    TD <- round(TD, 4)
    
    # Recurrent events
    if(recurrent == TRUE){
      TH <- sim_event(
        A=Ai, X=Xi, b1=b1i, b2=b2i, lambda=lambdaH0,
        k1=k1, k2=k2, alpha=gamma, beta1=beta1, beta2=beta2,
        recurrent = TRUE,start=start, end=end)
      TH <- round(as.numeric(TH), 4)
    }else{
      TH <- NULL
    }
    
    tp <- t[2:5]
    # Censoring distribution
    if(censor == "independent"){
      TC <- sample(tp, 1, prob = c(0.15,0.2,0.25, 0.4))
    }
    if(censor == "dependent"){
      Y1i <- Y1[i, 2:5]
      Y2i <- Y2[i, 2:5]
      probi <- expit( -1.0 + 0.8*Y1i - 0.5*Y2i )
      TCi <- rbinom(length(probi), 1, prob=probi)
      if(sum(TCi)==0){
        TC = tp[length(probi)]
      }else{
        TC = tp[min(which(TCi==1))]
      }
    }
    
    # Summarizing data without censoring
    TX <- min(tmax, TD)
    delta <- 1*(TD <= tmax)
    TH_before_death <- TH[TH < tmax]
    #TH_before_death <- TH[TH < TX]
    if(TD <= tmax){
      event_data_i <- data.frame(
        id = i,
        time = c(TH_before_death, TD, tmax),
        status = c( rep(2, length(TH_before_death)), 1, 0)
      )
    }else{
      event_data_i <- data.frame(
        id = i,
        time = c(TH_before_death, tmax),
        status = c( rep(2, length(TH_before_death)), 0)
      )
    }

    event_data <- rbind.data.frame(event_data,event_data_i)
    
    # Summarizing data with censoring
    TX <- min(TC, TD)
    delta <- 1*(TD <= TC)
    TH_before_death <- TH[TH < TC]
    if(TD <= TC){
      event_miss_i <- data.frame(
        id = i,
        time = c(TH_before_death, TD, TC),
        status = c( rep(2, length(TH_before_death)), 1, 0)
      )
    }else{
      event_miss_i <- data.frame(
        id = i,
        time = c(TH_before_death, TC),
        status = c( rep(2, length(TH_before_death)), 0)
      )
    }
    event_miss <- rbind.data.frame(event_miss,event_miss_i)
  }
  
  event_data <- event_data %>%
    distinct(id,time, status)
  
  event_miss <- event_miss %>%
    distinct(id, time, status)
  
  TX <- event_miss$time[event_miss$status ==0]
  #death <- event_miss$status[event_miss$status !=2]
  
  cov_miss <- cov_data %>%
    mutate(TX) %>%
    mutate(
      Y1_t3 = Y1_t3 * ifelse(TX < 3, NA, 1),
      Y1_t6 = Y1_t6 * ifelse(TX < 6, NA, 1),
      Y1_t9 = Y1_t9 * ifelse(TX < 9, NA, 1),
      Y1_t12 = Y1_t12 * ifelse(TX < 12, NA, 1),
      Y2_t3 = Y2_t3 * ifelse(TX < 3, NA, 1),
      Y2_t6 = Y2_t6 * ifelse(TX < 6, NA, 1),
      Y2_t9 = Y2_t9 * ifelse(TX < 9, NA, 1),
      Y2_t12 = Y2_t12 * ifelse(TX < 12, NA, 1)
    )
  
  return(list(
    cov_data=cov_data, event_data=event_data,
    cov_miss=cov_miss, event_miss=event_miss))
}
