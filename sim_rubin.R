#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Title:  Imputation of a mix of variable types                            #
# Description: This program contains the simulation for the imputation of  #
#              a mix of variable types                                     #    
# Author: Tuo Wang (email: twang437@wisc.edu, website: tuowang.rbind.io)   #
# Date:   Aug 2021                                                         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

#simid <- "10-20-ind"

# Preparation 

library(tidyverse) # data manipulation and data visualization
library(survival)  # fit survival models (e.g. Kaplan-meier)
library(reda)      # calculate non-parametric mean cumulative function
library(eha)       # fit parametric proportional hazard model
library(MASS)      # produce samples from multivariate normal distribution
library(glue)      # manipulate string
library(survRM2)

expit <- function(x){
  exp(x)/(1+exp(x))
}

# Load functions for data generation and imputation

source(file = here::here("Code",  "sim_data_fun.R"))  
source(file = here::here("Code",  "imputation_fun_new.R"))

#source(file="/ua/twang437/imputation-mix/sim/09_17/Code/sim_data_fun.R")
#source(file="/ua/twang437/imputation-mix/sim/09_17/Code/imputation_fun_new.R")

n         <- 500                      # Sample size (total)
tmax      <- 12                       # Maximum follow-up time 
tt        <- c(0, 3, 6, 9, 12)        # pre-specified time points 
beta1     <- c(0.0, 1, -0.5)          # fixed effects for continuous var
beta2     <- c(0.0, 1, -0.5)          # fixed effects for binary var
k1        <- 0.5                      # time effect for continuous var
k2        <- 0.15                     # time effect for binary var
lambdaD0  <- 0.08                     # baseline hazard for time-to-event
lambdaH0  <- 0.13                     # baseline hazard for recurrent events
alpha     <- c(-0.7, 0.5, 0.5, -0.5)  # coefficients for time-to-event     
gamma     <- c(-0.7, 0.5, 0.5, -0.5)  # coefficients for recurrent events   
censor    <- "dependent"            # type of censoring distribution

random.samp.d <- TRUE # random sampling when impute time-to-event
random.samp.r <- TRUE # random sampling when impute recurrent events
adjust <- TRUE        # adjustment for parameters

# The covariates that are used at each imputation steps. 
# d:time-to-event, r:recurrent event, c: cts var, b: binary var.
# Z2_t1 : means the number of recurrent events in the first time interval
# Y1_t3 : means the value of longitudinal data at time t3

var.names <- list(
  d1 = list(
    d = c( "X", "Z2_t1","Y1_t3","Y2_t3"),
    r = c( "X", "Z1_t1", "Z2_t1","Y1_t3","Y2_t3"),
    c = c( "X", "Z1_t1", "Z2_t1","Y1_t3","Y2_t3","Z2_t2"),
    b = c( "X", "Z1_t1", "Z2_t1","Y1_t3","Y2_t3","Z2_t2","Y1_t6")
  ),
  d2 = list(
    d = c( "X", "Z2_t1","Z2_t2","Y1_t6","Y2_t6"),
    r = c( "X", "Z1_t2", "Z2_t1","Z2_t2","Y1_t6","Y2_t6"),
    c = c( "X", "Z1_t2", "Z2_t1","Z2_t2","Y1_t6","Y2_t6","Z2_t3"),
    b = c( "X", "Z1_t2", "Z2_t1","Z2_t2","Y1_t6","Y2_t6","Z2_t3", "Y1_t9")
  ),
  d3 = list(
    d = c( "X", "Z2_t1","Z2_t2","Z2_t3","Y1_t9","Y2_t9"),
    r = c( "X",  "Z1_t3", "Z2_t1","Z2_t2","Z2_t3","Y1_t9","Y2_t9"),
    c = c( "X",  "Z1_t3", "Z2_t1","Z2_t2","Z2_t3","Y1_t9","Y2_t9","Z2_t4"),
    b = c( "X",  "Z1_t3", "Z2_t1","Z2_t2","Z2_t3","Y1_t9","Y2_t9","Z2_t4","Y1_t12")
  )
)




# Number of replicates
ns <- 100
# Number of imputation
nimp <- 50

# time sequence for estimating survival probabilities and MCF
#t.seq <- (0:365)/30.417 
t.seq <- seq(0, 12, length.out=366)
l <- length(t.seq)
# Results matrix

survest.imp.A0 <- NULL
survest.imp.A1 <- NULL
survest.dat.A0 <- NULL
survest.dat.A1 <- NULL
survest.mis.A0 <- NULL
survest.mis.A1 <- NULL

mcfest.imp.A0 <- NULL
mcfest.imp.A1 <- NULL
mcfest.dat.A0 <- NULL
mcfest.dat.A1 <- NULL
mcfest.mis.A0 <- NULL
mcfest.mis.A1 <- NULL

Y1.df <- NULL
Y2.df <- NULL

nimp_real <- NULL

results.imp <- NULL

rmt.dat.A0 <- NULL
rmt.dat.A1 <- NULL
rmt.imp.A0 <- NULL
rmt.imp.A1 <- NULL
rmt.wse.imp.A0 <- NULL
rmt.wse.imp.A1 <- NULL

for (i in 1:ns) {
  
  cat(glue("** Simulated Data {i}"), "\n")
  cat("Generating data \n")
  
  # Generating non-censored and censored data
  dat.try <- tryCatch(
    {  dat.obj <- sim_data(
      n = n, tmax = tmax, t=tt,
      k1=k1, k2=k2, beta1=beta1, beta2=beta2,
      alpha=alpha, gamma = gamma,
      lambdaD0 = lambdaD0, lambdaH0 = lambdaH0,
      recurrent = TRUE, censor = censor, start = 0, end = 5000
    )},
    error = function(e){e})
  if(inherits(dat.try, "error")){cat("error in simulating data \n");next}
  
  # non-censored data but maximum follow-up time is 12.
  cov.data <- dat.obj$cov_data
  event.data <- dat.obj$event_data
  
  # id for control and treatment groups
  id.A0 <- cov.data$id[cov.data$A==0]
  id.A1 <- cov.data$id[cov.data$A==1]
  n0 <- length(id.A0)
  n1 <- length(id.A1)
  
  # censored data
  cov.miss <- dat.obj$cov_miss
  event.miss <- dat.obj$event_miss
  
  # divide the data according treatment group 
  event.miss.A0 <- event.miss[event.miss$id %in% id.A0,]
  cov.miss.A0 <- cov.miss[cov.miss$id %in% id.A0, ]
  
  event.miss.A1 <- event.miss[event.miss$id %in% id.A1,]
  cov.miss.A1 <- cov.miss[cov.miss$id %in% id.A1, ]
  
  # Begin Multiple Imputation
  survest.imp.A0.it <- NULL
  survest.imp.A1.it <- NULL
  survse.imp.A0.it  <- NULL
  survse.imp.A1.it  <- NULL
  
  mcfest.imp.A0.it <- NULL
  mcfest.imp.A1.it <- NULL
  mcfse.imp.A0.it  <- NULL
  mcfse.imp.A1.it  <- NULL
  
  Y1.df.it <- NULL
  Y2.df.it <- NULL
  
  cts12.est.A0 <- NULL
  cts12.se.A0  <- NULL
  cts12.est.A1 <- NULL
  cts12.se.A1  <- NULL
  bin12.est.A0 <- NULL
  bin12.se.A0  <- NULL
  bin12.est.A1 <- NULL
  bin12.se.A1  <- NULL
  
  rmt.imp.A0.it <- NULL
  rmt.imp.A1.it <- NULL
  
  rmt.wse.imp.A0.it <- NULL
  rmt.wse.imp.A1.it <- NULL
  
  mimp.try <- tryCatch({
    for (it in c(1:nimp)) {
      cat(glue("** Multiple imputation {it}"), "\n")
      # Imputation: impute the data for control and treatment group seperately
      imp.try <- tryCatch(
        {imp.obj.A0 <- imputation(
          event.miss = event.miss.A0, cov.miss = cov.miss.A0,
          var.names = var.names, adjust=adjust,
          random.samp.d = random.samp.d, random.samp.r = random.samp.r
        )
        
        imp.obj.A1 <- imputation(
          event.miss = event.miss.A1, cov.miss = cov.miss.A1,
          var.names = var.names, adjust=adjust, 
          random.samp.d = random.samp.d, random.samp.r = random.samp.r
        )},
        error = function(e){e})
      if(inherits(imp.try, "error")){cat("error in imputation \n");break}
      
      event.imp.A0 <- imp.obj.A0$event.imp
      cov.imp.A0 <- imp.obj.A0$cov.imp
      
      event.imp.A1 <- imp.obj.A1$event.imp
      cov.imp.A1 <- imp.obj.A1$cov.imp

      event.imp <- bind_rows(event.imp.A0, event.imp.A1) %>% arrange(id)
      cov.imp <- bind_rows(cov.imp.A0, cov.imp.A1) %>% arrange(id)
      
      # Calculate the RMST and the SE estimator
      d.event.imp <- event.imp %>%
        filter(status != 2) %>%
        group_by(id) %>%
        slice(1) %>%
        ungroup()
      
      obj <- rmst2(d.event.imp$time, d.event.imp$status, cov.imp$A, tau=12)
      rmt.imp.A0.it    <- c(rmt.imp.A0.it, obj$RMST.arm0$rmst[1])
      rmt.imp.A1.it    <- c(rmt.imp.A1.it, obj$RMST.arm1$rmst[1])
      rmt.wse.imp.A0.it <- c(rmt.wse.imp.A0.it, obj$RMST.arm0$rmst[2])
      rmt.wse.imp.A1.it <- c(rmt.wse.imp.A1.it, obj$RMST.arm1$rmst[2])
      
      # Summarizing results!!
      # Death - imputed data
      d.event.imp.A0 <- event.imp %>%
        filter( id %in% id.A0 & status != 2) %>%
        group_by(id) %>%
        slice(1) %>%
        ungroup()
      
      d.event.imp.A1 <- event.imp %>%
        filter( id %in% id.A1 & status != 2) %>%
        group_by(id) %>%
        slice(1) %>%
        ungroup()
      
      km.imp.A0 <- survfit(Surv(time, status)~1,data=d.event.imp.A0)
      km.imp.A1 <- survfit(Surv(time, status)~1,data=d.event.imp.A1)
      
      survest.imp.A0.it <- rbind(survest.imp.A0.it, summary(km.imp.A0, times = t.seq)$surv)
      survest.imp.A1.it <- rbind(survest.imp.A1.it, summary(km.imp.A1, times = t.seq)$surv)
    
      # standard error estimator
      survse.stepfun.A0 <- stepfun(km.imp.A0$time, c(0,km.imp.A0$surv*km.imp.A0$std.err))
      survse.stepfun.A1 <- stepfun(km.imp.A1$time, c(0,km.imp.A1$surv*km.imp.A1$std.err))
      survse.imp.A0.it <- rbind(survse.imp.A0.it, survse.stepfun.A0(t.seq))
      survse.imp.A1.it <- rbind(survse.imp.A1.it, survse.stepfun.A1(t.seq))

      # Recurrent events - imputed data
      r.event.imp.A0 <- event.imp %>%
        add_row(id=1:n, time=12,status=0) %>%
        arrange(id) %>%
        distinct(id,time, status) %>%
        filter( id %in% id.A0 ) %>%
        filter( status != 1) %>%
        mutate( status = case_when(
          status == 0 ~ 0,
          status == 2 ~ 1
        ))
      
      r.event.imp.A1 <-   event.imp %>%
        add_row(id=1:n, time=12,status=0) %>%
        arrange(id) %>%
        distinct(id,time, status) %>%
        filter( id %in% id.A1 ) %>%
        filter( status != 1) %>%
        mutate( status = case_when(
          status == 0 ~ 0,
          status == 2 ~ 1
        ))
      
      mcf.imp.A0 <- mcf(Recur(time, id, status) ~ 1, data = r.event.imp.A0)
      mcf.imp.A1 <- mcf(Recur(time, id, status) ~ 1, data = r.event.imp.A1)
      
      mcf.stepfun.imp.A0 <- stepfun(mcf.imp.A0@MCF$time, c(0,mcf.imp.A0@MCF$MCF))
      mcf.stepfun.imp.A1 <- stepfun(mcf.imp.A1@MCF$time, c(0,mcf.imp.A1@MCF$MCF))
      
      mcf.imp.A0@MCF$MCF - 1.96 * mcf.imp.A0@MCF$se
      
      mcfest.imp.A0.it <- rbind(mcfest.imp.A0.it, mcf.stepfun.imp.A0(t.seq))
      mcfest.imp.A1.it <- rbind(mcfest.imp.A1.it, mcf.stepfun.imp.A1(t.seq))
      
      # standard error estimator
      mcfse.stepfun.A0 <- stepfun(mcf.imp.A0@MCF$time, c(0,mcf.imp.A0@MCF$se))
      mcfse.stepfun.A1 <- stepfun(mcf.imp.A1@MCF$time, c(0,mcf.imp.A1@MCF$se))
      mcfse.imp.A0.it <- rbind(mcfse.imp.A0.it, mcfse.stepfun.A0(t.seq))
      mcfse.imp.A1.it <- rbind(mcfse.imp.A1.it, mcfse.stepfun.A1(t.seq))
      
      # Continuous variable
      Y1_t12.A0 <- cov.imp$Y1_t12[cov.imp$A == 0]
      cts12.est.A0 <- c(cts12.est.A0, mean(Y1_t12.A0))
      cts12.se.A0  <- c(cts12.se.A0, sd(Y1_t12.A0)/sqrt(n0))
      
      Y1_t12.A1 <- cov.imp$Y1_t12[cov.imp$A == 1]
      cts12.est.A1 <- c(cts12.est.A1, mean(Y1_t12.A1))
      cts12.se.A1  <- c(cts12.se.A1, sd(Y1_t12.A1)/sqrt(n1))
      
      # Binary variable
      phat.A0 <- mean(cov.imp$Y2_t12[cov.imp$A == 0])
      bin12.est.A0 <- c(bin12.est.A0, phat.A0)
      bin12.se.A0  <- c(bin12.se.A0, sqrt(phat.A0*(1-phat.A0)/n0))
      
      phat.A1 <- mean(cov.imp$Y2_t12[cov.imp$A == 1])
      bin12.est.A1 <- c(bin12.est.A1, phat.A1)
      bin12.se.A1  <- c(bin12.se.A1, sqrt(phat.A1*(1-phat.A1)/n1))
      
      Y1.df.i <- cov.imp %>% 
        group_by(A) %>%
        summarise(
          across(starts_with("Y1"), ~mean(.x, na.rm=TRUE)),
          .groups = 'drop') %>%
        mutate(num.it = it) %>%
        pivot_longer(
          starts_with("Y"), names_prefix = "Y1_t",
          names_to = "time", values_to = "Y1") %>%
        mutate(type="Imputation") 
      
      Y1.df.it <- rbind.data.frame(Y1.df.it, Y1.df.i)
      
      Y2.df.i <- cov.imp %>% 
        group_by(A) %>%
        summarise(
          across(starts_with("Y2"), ~mean(.x, na.rm=TRUE)),
          .groups = 'drop') %>%
        mutate(num.it = it) %>%
        pivot_longer(
          starts_with("Y"), names_prefix = "Y2_t",
          names_to = "time", values_to = "Y2") %>%
        mutate(type="Imputation") 
      
      Y2.df.it <- rbind.data.frame(Y2.df.it, Y2.df.i)
    }
  },
  error=function(e){e})
  if(inherits(mimp.try, "error")){cat("error in imputation \n");next}
  
  if(is.null(survest.imp.A0.it)){next}
  
  # Calculate the RMST for non-censored data
  d.event.dat <- event.data %>%
    filter(status != 2) %>%
    group_by(id) %>%
    slice(1) %>%
    ungroup()
  
  obj <- rmst2(d.event.dat$time, d.event.dat$status, cov.data$A, tau=12)
  rmt.dat.A0 <- c(rmt.dat.A0, obj$RMST.arm0$rmst[1])
  rmt.dat.A1 <- c(rmt.dat.A1, obj$RMST.arm1$rmst[1])
  
  # censored data
  cov.miss <- dat.obj$cov_miss
  event.miss <- dat.obj$event_miss
  
  # divide the data according treatment group 
  event.miss.A0 <- event.miss[event.miss$id %in% id.A0,]
  cov.miss.A0 <- cov.miss[cov.miss$id %in% id.A0, ]
  
  event.miss.A1 <- event.miss[event.miss$id %in% id.A1,]
  cov.miss.A1 <- cov.miss[cov.miss$id %in% id.A1, ]
  
  # Death - non-censored data
  d.event.dat.A0 <- event.data %>%
    filter( id %in% id.A0 & status != 2) %>%
    group_by(id) %>%
    slice(1) %>%
    ungroup()
  
  d.event.dat.A1 <- event.data %>%
    filter( id %in% id.A1 & status != 2) %>%
    group_by(id) %>%
    slice(1) %>%
    ungroup()
  
  km.dat.A0 <- survfit(Surv(time, status)~1,data=d.event.dat.A0)
  km.dat.A1 <- survfit(Surv(time, status)~1,data=d.event.dat.A1)
  
  survest.dat.A0 <- rbind(survest.dat.A0, summary(km.dat.A0, times = t.seq)$surv)
  survest.dat.A1 <- rbind(survest.dat.A1, summary(km.dat.A1, times = t.seq)$surv)
  
  # Death - censored data
  d.event.mis.A0 <- event.miss %>%
    filter( id %in% id.A0 & status != 2) %>%
    group_by(id) %>%
    slice(1) %>%
    ungroup()
  
  d.event.mis.A1 <- event.miss %>%
    filter( id %in% id.A1 & status != 2) %>%
    group_by(id) %>%
    slice(1) %>%
    ungroup()
  
  km.mis.A0 <- survfit(Surv(time, status)~1,data=d.event.mis.A0)
  km.mis.A1 <- survfit(Surv(time, status)~1,data=d.event.mis.A1)
  
  survest.mis.A0 <- rbind(survest.mis.A0, summary(km.mis.A0, times = t.seq)$surv)
  survest.mis.A1 <- rbind(survest.mis.A1, summary(km.mis.A1, times = t.seq)$surv)
  
  
  # Recurrent events - non-censored data
  r.event.dat.A0 <- event.data %>%
    filter( id %in% id.A0 ) %>%
    filter( status != 1) %>%
    mutate( status = case_when(
      status == 0 ~ 0,
      status == 2 ~ 1
    ))
  r.event.dat.A1 <- event.data %>%
    filter( id %in% id.A1 ) %>%
    filter( status != 1) %>%
    mutate( status = case_when(
      status == 0 ~ 0,
      status == 2 ~ 1
    ))
  
  mcf.dat.A0 <- mcf(Recur(time, id, status) ~ 1, data = r.event.dat.A0)
  mcf.dat.A1 <- mcf(Recur(time, id, status) ~ 1, data = r.event.dat.A1)
  
  mcf.stepfun.dat.A0 <- stepfun(mcf.dat.A0@MCF$time, c(0,mcf.dat.A0@MCF$MCF))
  mcf.stepfun.dat.A1 <- stepfun(mcf.dat.A1@MCF$time, c(0,mcf.dat.A1@MCF$MCF))
  
  mcfest.dat.A0 <- rbind(mcfest.dat.A0, mcf.stepfun.dat.A0(t.seq))
  mcfest.dat.A1 <- rbind(mcfest.dat.A1, mcf.stepfun.dat.A1(t.seq))
  
  # Recurrent events - censored data
  r.event.mis.A0 <- event.miss %>%
    filter( id %in% id.A0 ) %>%
    filter( status != 1) %>%
    mutate( status = case_when(
      status == 0 ~ 0,
      status == 2 ~ 1
    ))
  
  r.event.mis.A1 <- event.miss %>%
    filter( id %in% id.A1 ) %>%
    filter( status != 1) %>%
    mutate( status = case_when(
      status == 0 ~ 0,
      status == 2 ~ 1
    ))
  
  mcf.mis.A0 <- mcf(Recur(time, id, status) ~ 1, data = r.event.mis.A0)
  mcf.mis.A1 <- mcf(Recur(time, id, status) ~ 1, data = r.event.mis.A1)
  
  mcf.stepfun.mis.A0 <- stepfun(mcf.mis.A0@MCF$time, c(0,mcf.mis.A0@MCF$MCF))
  mcf.stepfun.mis.A1 <- stepfun(mcf.mis.A1@MCF$time, c(0,mcf.mis.A1@MCF$MCF))
  
  mcfest.mis.A0 <- rbind(mcfest.mis.A0, mcf.stepfun.mis.A0(t.seq))
  mcfest.mis.A1 <- rbind(mcfest.mis.A1, mcf.stepfun.mis.A1(t.seq))
  
  # Restricted time
  
  rmt.imp.A0 <- c(rmt.imp.A0, mean(rmt.imp.A0.it))
  rmt.imp.A1 <- c(rmt.imp.A1, mean(rmt.imp.A1.it))
  
  rmt.wse.imp.A0 <- c(rmt.wse.imp.A0, sqrt(mean(rmt.wse.imp.A0.it^2)))
  rmt.wse.imp.A1 <- c(rmt.wse.imp.A1, sqrt(mean(rmt.wse.imp.A1.it^2)))
  
  # For time 12. 
  l <- length(t.seq)
  
  surv12.A0 <- data.frame(
    variable = "surv12",
    iteration = i,
    A = 0,
    est = mean(survest.imp.A0.it[,l]),
    sd_within = sqrt(mean(survse.imp.A0.it[,l]^2)),
    se_rubin = rubin_formula(survest.imp.A0.it[,l], survse.imp.A0.it[,l])
  )
  
  surv12.A1 <- data.frame(
    variable = "surv12",
    iteration = i,
    A = 1,
    est = mean(survest.imp.A1.it[,l]),
    sd_within = sqrt(mean(survse.imp.A1.it[,l]^2)),
    se_rubin = rubin_formula(survest.imp.A1.it[,l], survse.imp.A1.it[,l])
  )
  
  mcf12.A0 <- data.frame(
    variable = "mcf12",
    iteration = i,
    A = 0,
    est = mean(mcfest.imp.A0.it[,l]),
    sd_within = sqrt(mean(mcfse.imp.A0.it[,l]^2)),
    se_rubin = rubin_formula(mcfest.imp.A0.it[,l], mcfse.imp.A0.it[,l])
  )
  
  mcf12.A1 <- data.frame(
    variable = "mcf12",
    iteration = i,
    A = 1,
    est = mean(mcfest.imp.A1.it[,l]),
    sd_within = sqrt(mean(mcfse.imp.A1.it[,l]^2)),
    se_rubin = rubin_formula(mcfest.imp.A1.it[,l], mcfse.imp.A1.it[,l])
  )
  
  cts12.A0 <- data.frame(
    variable = "cts12",
    iteration = i,
    A = 0,
    est = mean(cts12.est.A0),
    sd_within = sqrt(mean(cts12.se.A0^2)),
    se_rubin = rubin_formula(cts12.est.A0, cts12.se.A0)
  )
  
  cts12.A1 <- data.frame(
    variable = "cts12",
    iteration = i,
    A = 1,
    est = mean(cts12.est.A1),
    sd_within = sqrt(mean(cts12.se.A1^2)),
    se_rubin = rubin_formula(cts12.est.A1, cts12.se.A1)
  )
  
  bin12.A0 <- data.frame(
    variable = "bin12",
    iteration = i,
    A = 0,
    est = mean(bin12.est.A0),
    sd_within = sqrt(mean(bin12.se.A0^2)),
    se_rubin = rubin_formula(bin12.est.A0, bin12.se.A0)
  )
  
  bin12.A1 <- data.frame(
    variable = "bin12",
    iteration = i,
    A = 1,
    est = mean(bin12.est.A1),
    sd_within = sqrt(mean(bin12.se.A1^2)),
    se_rubin = rubin_formula(bin12.est.A1, bin12.se.A1)
  )
  if(is.na(mean(bin12.est.A1))){stop("dd")}
  
  
  results.imp <- rbind.data.frame(
    results.imp,
    surv12.A0,
    surv12.A1,
    mcf12.A0,
    mcf12.A1,
    cts12.A0,
    cts12.A1,
    bin12.A0,
    bin12.A1
  )
  
  nimp_real <- c(nimp_real, nrow(survest.imp.A0.it))
  
  survest.imp.A0 <- rbind(survest.imp.A0, apply(survest.imp.A0.it, 2, mean))
  survest.imp.A1 <- rbind(survest.imp.A1, apply(survest.imp.A1.it, 2, mean))
  
  mcfest.imp.A0 <- rbind(mcfest.imp.A0, apply(mcfest.imp.A0.it, 2, mean))
  mcfest.imp.A1 <- rbind(mcfest.imp.A1, apply(mcfest.imp.A1.it, 2, mean))
  
  Y1.df.imp <- Y1.df.it %>%
    group_by(type, A, time) %>%
    summarise(
      Y1 = mean(Y1),
      .groups = 'drop'
    ) %>% 
    mutate(iteration = i) %>%
    dplyr::select(type, A, iteration, time ,Y1)
  
  Y2.df.imp <- Y2.df.it %>%
    group_by(type, A, time) %>%
    summarise(
      Y2 = mean(Y2),
      .groups = 'drop'
    ) %>% 
    mutate(iteration = i) %>%
    dplyr::select(type, A, iteration, time ,Y2)
  
  Y1.df.other <- cov.miss %>%
    mutate(type = "No imputation") %>%
    bind_rows(cov.data %>% mutate(type="No censoring")) %>%
    group_by(type, A) %>%
    summarise(
      across(starts_with("Y1"), ~mean(.x, na.rm=TRUE)),
      .groups = 'drop') %>%
    mutate(iteration = i) %>%
    pivot_longer(
      starts_with("Y"), names_prefix = "Y1_t",
      names_to = "time", values_to = "Y1")
  
  
  Y2.df.other <- cov.miss %>%
    mutate(type = "No imputation") %>%
    bind_rows(cov.data %>% mutate(type="No censoring")) %>%
    group_by(type, A) %>%
    summarise(
      across(starts_with("Y2"), ~mean(.x, na.rm=TRUE)),
      .groups = 'drop') %>%
    mutate(iteration = i) %>%
    pivot_longer(
      starts_with("Y"), names_prefix = "Y2_t",
      names_to = "time", values_to = "Y2")
  
  Y1.df <- rbind.data.frame(
    Y1.df, bind_rows(Y1.df.other, Y1.df.imp))
  
  Y2.df <- rbind.data.frame(
    Y2.df, bind_rows(Y2.df.other, Y2.df.imp))

  #save.image(file = glue("./sim{simid}.RData"))
  
}

mean(rmt.dat.A0)
mean(rmt.dat.A1)
mean(rmt.imp.A0)
mean(rmt.imp.A1)

sd(rmt.dat.A0)
sd(rmt.dat.A1)
mean(rmt.wse.imp.A0)
mean(rmt.wse.imp.A1)

#save.image(file = glue("./sim{simid}.RData"))

surv12.A0.TRUE <- mean(survest.dat.A0[,l])
surv12.A1.TRUE <- mean(survest.dat.A1[,l])
mcf12.A0.TRUE <- mean(mcfest.dat.A0[,l])
mcf12.A1.TRUE <- mean(mcfest.dat.A1[,l])
Y1_t12.A0.TRUE <- -0.00154
Y1_t12.A1.TRUE <- -0.498
Y2_t12.A0.TRUE <- 0.500
Y2_t12.A1.TRUE <- 0.353

cp <- function(est.true, est, se){
  1*( est.true >= est-1.96*se & est.true <= est+1.96*se )
}

results.imp$est.true <- rep(
  c(surv12.A0.TRUE,surv12.A1.TRUE,mcf12.A0.TRUE, mcf12.A1.TRUE,
    Y1_t12.A0.TRUE,Y1_t12.A1.TRUE,Y2_t12.A0.TRUE,Y2_t12.A1.TRUE),
  2000)

results.imp %>%
  mutate(cp = cp(est.true, est, se_rubin)) %>%
  group_by(variable, A) %>%
  summarise(
    mean(est,na.rm=TRUE), 
    sd(est,na.rm=TRUE), 
    mean(se_rubin,na.rm=TRUE), 
    mean(sd_within,na.rm=TRUE),mean(cp,na.rm = TRUE)
    ) %>% View()

cat(glue( "surv12 A0, EST:{round(mean(survest.dat.A0[,l]),6)}, SD:{round(sd(survest.dat.A0[,l]),6)}
          surv12 A1, EST:{round(mean(survest.dat.A1[,l]),6)}, SD:{round(sd(survest.dat.A1[,l]),6)},
          mcf12 A0, EST:{round(mean(mcfest.dat.A0[,l]),6)}, SD:{round(sd(mcfest.dat.A0[,l]),6)},
          mcf12 A1, EST:{round(mean(mcfest.dat.A1[,l]),6)}, SD:{round(sd(mcfest.dat.A1[,l]),6)}" ))

Y1.df %>%
  filter(time == 12) %>%
  group_by(type, A) %>%
  summarise(mean(Y1), sd(Y1))

Y2.df %>%
  filter(time == 12) %>%
  group_by(type, A) %>%
  summarise(mean(Y2), sd(Y2))

Y1.df %>%
  mutate(time = as.numeric(time)) %>%
  mutate(A = case_when(
    A == 0 ~ "Control",
    A == 1 ~ "Treatment"
  )) %>%
  group_by(type, A, time) %>%
  summarise(
    Y1_EST = mean(Y1),
    Y1_SD = sd(Y1),
    .groups = 'drop'
  ) %>%
  filter(time == 12)

Y2.df %>%
  mutate(time = as.numeric(time)) %>%
  mutate(A = case_when(
    A == 0 ~ "Control",
    A == 1 ~ "Treatment"
  )) %>%
  group_by(type, A, time) %>%
  summarise(
    Y2_EST = mean(Y2),
    Y2_SD = sd(Y2),
    .groups = 'drop'
  ) %>%
  filter(time == 12)

#load(file=here::here("Server/10_18b/sim2-3.RData"))

dim(survest.imp.A0)
avg.survest.imp.A0 <- apply(survest.imp.A0, 2, mean)
avg.survest.imp.A1 <- apply(survest.imp.A1, 2, mean)

avg.survest.dat.A0 <- apply(survest.dat.A0, 2, mean)
avg.survest.dat.A1 <- apply(survest.dat.A1, 2, mean)

avg.survest.mis.A0 <- apply(survest.mis.A0, 2, mean)
avg.survest.mis.A1 <- apply(survest.mis.A1, 2, mean)

par(mfrow=c(1,2))
plot( t.seq, avg.survest.imp.A0,col=2, type='l',xlim=c(0,12), ylim =c(0,1),lwd=2,lty=1)
lines(t.seq, avg.survest.dat.A0,col=3, type='l',xlim=c(0,12), ylim =c(0,1),lwd=2,lty=1)
lines(t.seq, avg.survest.mis.A0,col=4, type='l',xlim=c(0,12), ylim =c(0,1),lwd=2,lty=1)

plot( t.seq, avg.survest.imp.A1,col=2, type='l',xlim=c(0,12), ylim =c(0,1),lwd=2,lty=1)
lines(t.seq, avg.survest.dat.A1,col=3, type='l',xlim=c(0,12), ylim =c(0,1),lwd=2,lty=1)
lines(t.seq, avg.survest.mis.A1,col=4, type='l',xlim=c(0,12), ylim =c(0,1),lwd=2,lty=1)

# Plot MCF


avg.mcfest.imp.A0 <- apply(mcfest.imp.A0, 2, mean)
avg.mcfest.imp.A1 <- apply(mcfest.imp.A1, 2, mean)

avg.mcfest.dat.A0 <- apply(mcfest.dat.A0, 2, mean)
avg.mcfest.dat.A1 <- apply(mcfest.dat.A1, 2, mean)

avg.mcfest.mis.A0 <- apply(mcfest.mis.A0, 2, mean)
avg.mcfest.mis.A1 <- apply(mcfest.mis.A1, 2, mean)

par(mfrow=c(1,2))
plot( t.seq, avg.mcfest.imp.A0,col=2, type='l',xlim=c(0,12),lwd=2,lty=1)
lines(t.seq, avg.mcfest.dat.A0,col=3, type='l',xlim=c(0,12),lwd=2,lty=1)
lines(t.seq, avg.mcfest.mis.A0,col=4, type='l',xlim=c(0,12),lwd=2,lty=1)

plot( t.seq, avg.mcfest.imp.A1,col=2, type='l',xlim=c(0,12),lwd=2,lty=1)
lines(t.seq, avg.mcfest.dat.A1,col=3, type='l',xlim=c(0,12),lwd=2,lty=1)
lines(t.seq, avg.mcfest.mis.A1,col=4, type='l',xlim=c(0,12),lwd=2,lty=1)

# Plot continuous variable

col_palette <-  c("#386cb0","#fdb462","#7fc97f")
Y1.df %>%
  mutate(time = as.numeric(time)) %>%
  mutate(A = case_when(
    A == 0 ~ "Control",
    A == 1 ~ "Treatment"
  )) %>%
  group_by(type, A, time) %>%
  summarise(
    Y1 = mean(Y1),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x=time, y=Y1, color=type )) +
  geom_point(aes(shape=type),size=2.5,alpha=0.7) +
  geom_line(aes(linetype=type),size=1.0, alpha=0.7) +
  scale_color_manual(values = col_palette) +
  scale_linetype_manual(values = c(2,3,4)) +
  scale_x_continuous(breaks = c(0,3,6,9,12)) +
  facet_wrap(~A) +
  theme_bw()

# Plot binary variable

Y2.df %>%
  mutate(time = as.numeric(time)) %>%
  mutate(A = case_when(
    A == 0 ~ "Control",
    A == 1 ~ "Treatment"
  )) %>%
  group_by(type, A, time) %>%
  summarise(
    Y2 = mean(Y2),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x=time, y=Y2, color=type )) +
  geom_point(aes(shape=type),size=2.5,alpha=0.7) +
  geom_line(aes(linetype=type),size=1.0, alpha=0.7) +
  scale_color_manual(values = col_palette) +
  scale_linetype_manual(values = c(2,3,4)) +
  scale_x_continuous(breaks = c(0,3,6,9,12)) +
  facet_wrap(~A) +
  theme_bw()
 
