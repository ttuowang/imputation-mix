#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Title:  Imputation of a mix of variable types                            #
# Description: This program contains the imputation program of             #
#              a mix of variable types                                     #    
# Author: Tuo Wang (email: twang437@wisc.edu, website: tuowang.rbind.io)   #
# Date:   Jun, July & Aug 2021                                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 


rubin_formula <- function(est, se){
  m <- length(est)
  #est.avg <- mean(est)
  var.w <- mean(se^2)
  var.b <- var(est)
  var.t <- var.w + var.b + var.b/m
  se.t <- sqrt(var.t)
  return(se.t)
}



update_Z_mat <- function(event.miss, t){
  
  #' event.miss: the long format of the event history data. status=0: censored, 
  #'             status=1: time-to-event, status=2: recurrent event
  #' t: pre-specified time-points
  #' return a matrix contains all covariates being used in the imputation algorithm
  
  Z.mat <- event.miss %>%
    group_by(id) %>%
    mutate(death = sum((status==1)*1)) %>%
    mutate(
      TX = ifelse(sum(death)==0,tail(time,1),time[status==1]),
      TC = tail(time,1)
      ) %>% 
    mutate( 
      Z1_t1 = case_when(
        death == 0 & TX == 0 ~ NA_real_,
        death == 1 & TX <= t[1] ~ 1,
        TRUE ~ 0
      ),
      Z1_t2 = case_when(
        death == 0 & TX <= t[1] ~ NA_real_,
        death == 1 & TX <= t[2] ~ 1,
        TRUE ~ 0
      ),
      Z1_t3 = case_when(
        death == 0 & TX <= t[2] ~ NA_real_,
        death == 1 & TX <= t[3] ~ 1,
        TRUE ~ 0
      ),
      Z1_t4 = case_when(
        death == 0 & TX <= t[3] ~ NA_real_,
        death == 1 & TX <= t[4] ~ 1,
        TRUE ~ 0
      )
    ) %>%
    mutate(
      Z2_t1 = case_when(
        TC == 0 ~ NA_real_,
        TRUE ~ 1*(status==2 & time>=0    & time < t[1])
      ),
      Z2_t2 = case_when(
        TC <= t[1] ~ NA_real_,
        TRUE ~ 1*(status==2 & time>=t[1] & time < t[2])
      ),
      Z2_t3 = case_when(
        TC <= t[2] ~ NA_real_,
        TRUE ~ 1*(status==2 & time>=t[2] & time < t[3])
      ),
      Z2_t4 = case_when(
        TC <= t[3] ~ NA_real_,
        TRUE ~ 1*(status==2 & time>=t[3] & time <= t[4])
      )
    ) %>% 
    mutate(
      Z2_t1 = sum(Z2_t1),
      Z2_t2 = sum(Z2_t2),
      Z2_t3 = sum(Z2_t3),
      Z2_t4 = sum(Z2_t4)
    ) %>%
    slice(1) %>%
    ungroup() #%>%
    # filter(status != 2)
  
  return(Z.mat)
}

update_Z_lst <- function(cov.miss, Z.mat){
  
  # Given the covariates data and the Z matrix, return the Z^(j) matrix for each
  # time points
  
  Z.t1 <- cov.miss %>%
    dplyr::select(id, A, X, Y1_t0, Y2_t0, Y1_t3, Y2_t3) %>%
    mutate(Z1_t1 = Z.mat$Z1_t1, Z2_t1 = Z.mat$Z2_t1) %>%
    relocate(id, A, X, Y1_t0, Y2_t0, Z1_t1, Z2_t1)
  
  Z.t2 <- cov.miss %>%
    dplyr::select(id, A, X, Y1_t0, Y2_t0, Y1_t6, Y2_t6) %>%
    mutate(Z1_t2 = Z.mat$Z1_t2, Z2_t2 = Z.mat$Z2_t2) %>%
    relocate(id, A, X, Y1_t0, Y2_t0, Z1_t2, Z2_t2)
  
  Z.t3 <- cov.miss %>%
    dplyr::select(id, A, X, Y1_t0, Y2_t0, Y1_t9, Y2_t9) %>%
    mutate(Z1_t3 = Z.mat$Z1_t3, Z2_t3 = Z.mat$Z2_t3) %>%
    relocate(id, A, X, Y1_t0, Y2_t0, Z1_t3, Z2_t3)
  
  Z.t4 <- cov.miss %>%
    dplyr::select(id, A, X, Y1_t0, Y2_t0, Y1_t12, Y2_t12) %>%
    mutate(Z1_t4 = Z.mat$Z1_t4, Z2_t4 = Z.mat$Z2_t4) %>%
    relocate(id, A, X, Y1_t0, Y2_t0, Z1_t4, Z2_t4)
  
  Z <- list( t1 = Z.t1,
             t2 = Z.t2,
             t3 = Z.t3,
             t4 = Z.t4 )
  return(Z)
}

imputation <- function(
  event.miss, cov.miss, 
  var.names,
  adjust=TRUE, 
  random.samp.d = FALSE, 
  random.samp.r = FALSE,
  verbose=FALSE, 
  verbose2=TRUE){
  
  if(verbose2){
    cat("begin imputation \n")
  }
  
  t <- c(3,6,9,12) # Fixed time points
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ##   Begin imputation for (3, 6]           ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  k <- 1
  tk <- t[k]
  
  if(verbose2){
    cat(glue("Imputation for ({tk}, {tk+3}]"),"\n")
  }
  
  # impute time-to-event time
  ## Update Z matrix, death and censoring time

  Z.mat <- update_Z_mat(event.miss, t)
  Z.lst <- update_Z_lst(cov.miss, Z.mat)
  Z.t1 <- Z.lst$t1
  Z.t2 <- Z.lst$t2
  Z.t3 <- Z.lst$t3
  Z.t4 <- Z.lst$t4
  
  TX <- Z.mat$TX
  death <- Z.mat$death
  
  Z.df <- Z.t1 %>%
    left_join(
      Z.t2 %>%
        dplyr::select(id, Z1_t2, Z2_t2, Y1_t6, Y2_t6),
      by = "id"
    ) %>%
    mutate(TX, death) %>%
    filter( !(TX <=tk & death==1))
  
  target.var <- glue("Z1_t{k+1}")
  
  # Z.df.obs <- Z.df %>%
  #   filter(!is.na(!!sym(target.var))) %>%
  #   mutate(time = TX-tk,
  #          status = death)
  
  Z.df.obs <- Z.df %>%
    filter(!is.na(!!sym(target.var))) %>%
    mutate(time = pmin(TX - tk, 3),
           status = ifelse(TX - tk < 3, 1, 0))
  
  Z.df.mis <- Z.df %>%
    filter(is.na(!!sym(target.var)))
  
  n.obs  <- nrow(Z.df.obs)
  id.obs <- Z.df.obs$id
  n.mis  <- nrow(Z.df.mis)
  id.mis <- Z.df.mis$id
  
  formula1 <- as.formula(
    glue(
      "Surv(time, status)", " ~ ",
      glue_collapse(var.names[[k]]$d,sep = " + ")
    )
  )
  
  
  fit1 <- survreg(
    formula1,
    dist='exponential',
    data = Z.df.obs)
  
  # Fit a parametric proportional model for imputing time-to-event
  # fit1 <- phreg(
  #   formula1, 
  #   dist = "weibull", shape = 1,
  #   data = Z.df.obs)
  
  theta1.mu <- fit1$coefficients
  theta1.var <- fit1$var
  if(random.samp.d){
    theta1.tilde <- MASS::mvrnorm(1, theta1.mu, theta1.var)
  }else{
    theta1.tilde <- theta1.mu 
  }
  
  event.miss.tmp <- event.miss
  for (i in c(1:n.mis)) {
    idi <- id.mis[i]
    coef.tilde <- theta1.tilde
    cov.mati <- Z.df.mis %>%
      filter(id == idi) %>%
      dplyr::select(all_of(var.names[[k]]$d)) %>%
      as.matrix()
    #cov.mati <- cbind(cov.mati,-1)
    cov.mati <- cbind(-1, -cov.mati)
    ratei <- exp(
      cov.mati %*% coef.tilde - 
        adjust * cov.mati %*% theta1.var %*% t(cov.mati)/2)

    Y1.imp <- rexp(1, rate = ratei)
    Z1.imp <- 1*(Y1.imp <= 3)
    #Z1.imp <- ifelse(is.na(Z1.imp), 0, Z1.imp)
    Z.t2 <- Z.t2 %>%
      rows_update(tibble(id=idi, Z1_t2=Z1.imp), by="id")
    if(!is.na(Z1.imp) & Z1.imp == 1){
      event.miss.tmp <- event.miss.tmp %>%
        rows_delete(
          tibble(id=idi,time=tk, status=0),
          by=c("id","time","status")) %>%
        add_row(id=idi, time=Y1.imp+tk, status=1) %>%
        add_row(id=idi, time=tk+3, status=0)
      if(verbose){
        cat(glue(
          "Impute {round(Y1.imp+tk,4)} death time for subject {idi} ",
          "at ({tk},{tk+3}]"), "\n")
      }
    }else{
      event.miss.tmp <- event.miss.tmp %>%
        rows_delete(
          tibble(id=idi,time=tk, status=0),
          by=c("id","time","status")) %>%
        add_row(id=idi, time=tk+3, status=0)
      if(verbose){
        cat(glue(
          "Impute {tk+3} censor time for subject {idi} ",
          "at ({tk},{tk+3}]"), "\n")
      }
      
    }
  }
  event.miss.tmp <- event.miss.tmp %>%
    arrange(id, time)
  
  idzz <- event.miss.tmp$id[event.miss.tmp$time==tk & event.miss.tmp$status==0]
  
  event.miss.tmp <- event.miss.tmp %>%
    rows_delete(tibble(id=idzz,time=tk, status=0),
                by=c("id","time","status")) %>%
    add_row(id=idzz,time=tk+3, status=0) %>%
    arrange(id, time)
  
  ## Impute recurrent events time
  ## Update death and censoring time
  Z.mat <- update_Z_mat(event.miss.tmp, t)
  TX <- Z.mat$TX
  death <- Z.mat$death
  
  Z.df <- Z.t1 %>%
    left_join(
      Z.t2 %>%
        dplyr::select(id, Z1_t2, Z2_t2, Y1_t6, Y2_t6),
      by = "id"
    ) %>%
    mutate(TX, death) 
  
  target.var <- glue("Z2_t{k+1}")
  
  Z.df.obs <- Z.df %>%
    filter(!is.na(!!sym(target.var))) 
  
  Z.df.mis <- Z.df %>%
    filter(is.na(!!sym(target.var)))
  
  n.obs  <- nrow(Z.df.obs)
  n.mis  <- nrow(Z.df.mis)
  id.mis <- Z.df.mis$id
  
  formula2 <- as.formula(
    glue(
      glue("Z2","t{k+1}",.sep = "_"), " ~ ",
      glue_collapse(var.names[[k]]$r,sep = " + ")
    )
  )
  
  # Fit a negative binomial model for imputing recurrent events
  
  fit2 <- glm.nb(formula2, data = Z.df.obs)
  theta2.mu <- coef(fit2)
  theta2.var <- vcov(fit2)

  if(random.samp.r){
    theta2.tilde <- MASS::mvrnorm(1, theta2.mu, theta2.var)
  }else{
    theta2.tilde <- theta2.mu
  }
  
  event.miss.end <- event.miss.tmp
  inf.count <- 0
  for (i in c(1:n.mis)) {
    idi <- id.mis[i]
    TXi <- Z.df.mis$TX[i]
    deathi <- Z.df.mis$death[i]
    # 
    # Ai <- Z.df.mis$A[i]
    # Xi <- Z.df.mis$X[i]

    cov.mati <- Z.df.mis %>%
      filter(id == idi) %>%
      dplyr::select(all_of(var.names[[k]]$r)) %>%
      as.matrix()
    cov.mati <- cbind(1, cov.mati)
    lambda.t2 <- exp(
      cov.mati %*% theta2.tilde - 
        adjust * cov.mati %*% theta2.var %*% t(cov.mati)/2)/3
    if(lambda.t2 > 100){
      inf.count = inf.count +1}
    #lambda.t2 <- min(lambda.t2, 0.5) # Set a threshold 
    Yend <- tk
    dti <- ifelse(deathi == 1, TXi, tk+3)
    dti <- tk+3
    Z2.imp <- 0
    while (Yend <= min(tk+3,dti)) {
      Y2.imp <- rexp(1, lambda.t2)
      Yend <- Yend + Y2.imp
      if(Yend >= min(tk+3,dti)){break}
      event.miss.end <- event.miss.end %>%
        add_row(id=idi, time=Yend, status=2)
      Z2.imp <- Z2.imp + 1
    }
    if(verbose==TRUE){
      cat(glue(
        "Impute {Z2.imp} recurrent events for subject {idi} ",
        "with rate {round(lambda.t2,4)} ",
        "at ({tk},{tk+3}]"), "\n")
    }
    Z.t2 <- Z.t2 %>%
      rows_update(tibble(id=idi, Z2_t2=Z2.imp), by="id")
  }
  
  event.miss.end <- event.miss.end %>%
    arrange(id, time)
  
  ## Impute continuous variable
  Z.df <- Z.t1 %>%
    left_join(
      Z.t2 %>%
        dplyr::select(id, Z1_t2, Z2_t2, Y1_t6, Y2_t6),
      by = "id"
    ) %>%
    mutate(TX, death) 
  
  target.var <- glue("Y1_t{tk+3}")
  
  Z.df.obs <- Z.df %>%
    filter(!is.na(!!sym(target.var))) 
  
  Z.df.mis <- Z.df %>%
    filter(is.na(!!sym(target.var)))
  
  n.obs  <- nrow(Z.df.obs)
  n.mis  <- nrow(Z.df.mis)
  id.obs <- Z.df.obs$id
  id.mis <- Z.df.mis$id
  
  # Fit a linear regression model for continuous var
  
  formula3 <- as.formula(
    glue(
      "Y1_t{tk+3}", " ~ ",
      glue_collapse(var.names[[k]]$c,sep = " + ")
    )
  )
  
  fit3 <- lm(formula3, data = Z.df.obs, x=TRUE, y=TRUE)
  theta3.mu <- coef(fit3)
  theta3.var <- vcov(fit3)
  x <- fit3$x
  y <- fit3$y
  nj <- nrow(x)
  pj <- ncol(x)
  L <- chol(solve(t(x) %*% x))
  epsilon.sd <- sqrt(deviance(fit3)/df.residual(fit3))

  Zi <- MASS::mvrnorm(1, rep(0,pj), diag(1,nrow=pj))
  gi <- rchisq(1, df=nj-pj-1)
  sigma2i <- epsilon.sd^2*(nj-pj-1)/gi
  sigmai <- sqrt(sigma2i)
  theta3.tilde <- theta3.mu + sigmai*t(L)%*%Zi
  
  cov.mati <- Z.df.mis %>%
    dplyr::select(all_of(var.names[[k]]$c)) %>%
    as.matrix()
  cov.mati <- cbind(1, cov.mati)
  
  Y1_t6.imp <- cov.mati %*% theta3.tilde + rnorm(n.mis)*sigmai
  Z.t2$Y1_t6[Z.t2$id %in% id.mis] <- Y1_t6.imp
  cov.miss$Y1_t6[cov.miss$id %in% id.mis] <- Y1_t6.imp
  
  # Y1_t6.imp <- NULL
  # for (i in c(1:n.mis)) {
  #   idi <- id.mis[i]
  #   Zi <- MASS::mvrnorm(1, rep(0,pj), diag(1,nrow=pj))
  #   gi <- rchisq(1, df=nj-pj-1)
  #   sigma2i <- epsilon.sd^2*(nj-pj-1)/gi
  #   sigmai <- sqrt(sigma2i)
  #   theta3.tilde <- theta3.mu + sigmai*t(L)%*%Zi
  #   
  #   cov.mati <- Z.df.mis %>%
  #     filter(id == idi) %>%
  #     dplyr::select(all_of(var.names[[k]]$c)) %>%
  #     as.matrix()
  #   cov.mati <- cbind(1, cov.mati)
  #   Y1_t6.imp <- c(Y1_t6.imp, cov.mati %*% theta3.tilde + rnorm(1)*sigmai)
  # }
  # Z.t2$Y1_t6[Z.t2$id %in% id.mis] <- Y1_t6.imp
  # cov.miss$Y1_t6[cov.miss$id %in% id.mis] <- Y1_t6.imp
  # 
  # fit3 <- lm(formula3, data = Z.df.obs)
  # theta3.mu <- coef(fit3)
  # theta3.var <- vcov(fit3)
  # theta3.tilde <-  MASS::mvrnorm(n.mis, theta3.mu, theta3.var)
  # cov.mati <- Z.df.mis %>%
  #   dplyr::select(all_of(var.names[[k]]$c)) %>%
  #   as.matrix()
  # cov.mati <- cbind(1, cov.mati)
  # epsilon.sd <- sqrt(deviance(fit3)/df.residual(fit3))
  # epsilon.Y <- rnorm(n.mis, mean = 0, sd=epsilon.sd)
  # Y1_t6.imp <- apply(cov.mati * theta3.tilde,1,sum) + epsilon.Y
  # Z.t2$Y1_t6[Z.t2$id %in% id.mis] <- Y1_t6.imp
  # cov.miss$Y1_t6[cov.miss$id %in% id.mis] <- Y1_t6.imp
  
  # Impute binary variables
  Z.df <- Z.t1 %>%
    left_join(
      Z.t2 %>%
        dplyr::select(id, Z1_t2, Z2_t2, Y1_t6, Y2_t6),
      by = "id"
    ) %>%
    mutate(TX, death) 
  
  target.var <- glue("Y2_t{tk+3}")
  
  Z.df.obs <- Z.df %>%
    filter(!is.na(!!sym(target.var))) 
  
  Z.df.mis <- Z.df %>%
    filter(is.na(!!sym(target.var)))
  
  n.obs  <- nrow(Z.df.obs)
  n.mis  <- nrow(Z.df.mis)
  id.obs <- Z.df.obs$id
  id.mis <- Z.df.mis$id
  
  formula4 <- as.formula(
    glue(
      "Y2_t{tk+3}", " ~ ",
      glue_collapse(var.names[[k]]$b,sep = " + ")
    )
  )
  
  # Fit a logistic regression for imputing binary var
  
  fit4 <- glm(formula4, data = Z.df.obs, family = "binomial")
  theta4.mu <- coef(fit4)
  theta4.var <- vcov(fit4)
  
  theta4.tilde <-  MASS::mvrnorm(1, theta4.mu, theta4.var)
  cov.mati <- Z.df.mis %>%
    dplyr::select(all_of(var.names[[k]]$b)) %>%
    as.matrix()
  cov.mati <- cbind(1, cov.mati)
  Y2_t6.imp <- expit(cov.mati %*% theta4.tilde)# + epsilon.Y
  U <- runif(n.mis)
  #Y2_t6.imp <- (Y2_t6.imp>U)*1
  Y2_t6.imp <- (Y2_t6.imp > 0.5)*1
  cov.miss$Y2_t6[cov.miss$id %in% id.mis] <- Y2_t6.imp
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ### Next interval (6, 9] Impute death time ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  event.miss <- event.miss.end
  k <- 2
  tk <- t[k]
  
  if(verbose2){cat(glue("Imputation for ({tk}, {tk+3}]"),"\n")}
  
  # Update Z matrix, death time and censor status
  Z.mat <- update_Z_mat(event.miss, t)
  Z.lst <- update_Z_lst(cov.miss, Z.mat)
  Z.t1 <- Z.lst$t1
  Z.t2 <- Z.lst$t2
  Z.t3 <- Z.lst$t3
  Z.t4 <- Z.lst$t4
  
  TX <- Z.mat$TX
  death <- Z.mat$death
  
  Z.df <- Z.t1 %>%
    left_join(
      Z.t2 %>%
        dplyr::select(id, Z1_t2, Z2_t2, Y1_t6, Y2_t6),
      by = "id"
    ) %>%
    left_join(
      Z.t3 %>%
        dplyr::select(id, Z1_t3, Z2_t3, Y1_t9, Y2_t9),
      by = "id"
    ) %>%
    mutate(TX, death) %>%
    filter( !(TX <=tk & death==1))
  
  target.var <- glue("Z1_t{k+1}")
  
  # Z.df.obs <- Z.df %>%
  #   filter(!is.na(!!sym(target.var))) %>%
  #   mutate(time = TX - tk,
  #          status = death)
  
  Z.df.obs <- Z.df %>%
    filter(!is.na(!!sym(target.var)))%>%
    mutate(time = pmin(TX - tk, 3),
           status = ifelse(TX - tk < 3, 1, 0))
  
  Z.df.mis <- Z.df %>%
    filter(is.na(!!sym(target.var)))
  
  n.obs <- nrow(Z.df.obs)
  id.obs <- Z.df.obs$id
  
  n.mis <- nrow(Z.df.mis)
  id.mis <- Z.df.mis$id
  
  formula1 <- as.formula(
    glue(
      "Surv(time, status)", " ~ ",
      glue_collapse(var.names[[k]]$d,sep = " + ")
    )
  )
  
  fit1 <- survreg(
    formula1,
    dist='exponential',
    data = Z.df.obs)
  
  # fit1 <- phreg(
  #   formula1, 
  #   dist = "weibull", shape = 1,
  #   data = Z.df.obs)
  
  theta1.mu <- fit1$coefficients
  theta1.var <- fit1$var
  if(random.samp.d){
    theta1.tilde <- MASS::mvrnorm(1, theta1.mu, theta1.var)
  }else{
    theta1.tilde <- theta1.mu 
  }
  event.miss.tmp <- event.miss
  for (i in c(1:n.mis)) {
    idi <- id.mis[i]
    coef.tilde <- theta1.tilde
    cov.mati <- Z.df.mis %>%
      filter(id == idi) %>%
      dplyr::select(all_of(var.names[[k]]$d)) %>%
      as.matrix()
    #cov.mati <- cbind(cov.mati,-1)
    cov.mati <- cbind(-1, -cov.mati)
    ratei <- exp(
      cov.mati %*% coef.tilde - 
        adjust * cov.mati %*% theta1.var %*% t(cov.mati)/2)
    
    Y1.imp <- rexp(1, rate = ratei)
    Z1.imp <- 1*(Y1.imp <= 3)
    Z.t3 <- Z.t3 %>%
      rows_update(tibble(id=idi, Z1_t3=Z1.imp), by="id")
    if(!is.na(Z1.imp) & Z1.imp == 1){
      event.miss.tmp <- event.miss.tmp %>%
        rows_delete(
          tibble(id=idi,time=tk, status=0),
          by=c("id","time","status")) %>%
        add_row(id=idi, time=Y1.imp+tk, status=1) %>%
        add_row(id=idi, time=tk+3, status=0)
      if(verbose == TRUE){
        cat(glue(
          "Impute {round(Y1.imp+tk,4)} death time for subject {idi} ",
          "at ({tk},{tk+3}]"), "\n")
      }
    }else{
      event.miss.tmp <- event.miss.tmp %>%
        rows_delete(
          tibble(id=idi,time=tk, status=0),
          by=c("id","time","status")) %>%
        add_row(id=idi, time=tk+3, status=0)
      if(verbose == TRUE){
        cat(glue(
          "Impute {tk+3} censor time for subject {idi} ",
          "at ({tk},{tk+3}]"), "\n")
      }
    }
  }
  event.miss.tmp <- event.miss.tmp %>%
    arrange(id, time)
  
  idzz <- event.miss.tmp$id[event.miss.tmp$time==tk & event.miss.tmp$status==0]
  
  event.miss.tmp <- event.miss.tmp %>%
    rows_delete(tibble(id=idzz,time=tk, status=0),
                by=c("id","time","status")) %>%
    add_row(id=idzz,time=tk+3, status=0) %>%
    arrange(id, time)
  
  ## Impute recurrent time
  ## Update death and censoring time
  Z.mat <- update_Z_mat(event.miss.tmp, t)
  TX <- Z.mat$TX
  death <- Z.mat$death
  
  Z.df <- Z.t1 %>%
    left_join(
      Z.t2 %>%
        dplyr::select(id, Z1_t2, Z2_t2, Y1_t6, Y2_t6),
      by = "id"
    ) %>%
    left_join(
      Z.t3 %>%
        dplyr::select(id, Z1_t3, Z2_t3, Y1_t9, Y2_t9),
      by = "id"
    ) %>%
    mutate(TX, death) 
  
  target.var <- glue("Z2_t{k+1}")
  
  Z.df.obs <- Z.df %>%
    filter(!is.na(!!sym(target.var))) 
  
  Z.df.mis <- Z.df %>%
    filter(is.na(!!sym(target.var)))
  
  n.obs <- nrow(Z.df.obs)
  n.mis <- nrow(Z.df.mis)
  id.mis <- Z.df.mis$id
  
  formula2 <- as.formula(
    glue(
      glue("Z2","t{k+1}",.sep = "_"), " ~ ",
      glue_collapse(var.names[[k]]$r,sep = " + ")
    )
  )
  
  fit2 <- glm.nb(formula2, data = Z.df.obs)
  theta2.mu <- coef(fit2)
  theta2.var <- vcov(fit2)
  if(random.samp.r){
    theta2.tilde <- MASS::mvrnorm(1, theta2.mu, theta2.var)
  }else{
    theta2.tilde <- theta2.mu
  }
  event.miss.end <- event.miss.tmp
  for (i in c(1:n.mis)) {
    idi <- id.mis[i]
    TXi <- Z.df.mis$TX[i]
    deathi <- Z.df.mis$death[i]
    Ai <- Z.df.mis$A[i]
    
    cov.mati <- Z.df.mis %>%
      filter(id == idi) %>%
      dplyr::select(all_of(var.names[[k]]$r)) %>%
      as.matrix()
    cov.mati <- cbind(1, cov.mati)
    lambda.t2 <- exp(
      cov.mati %*% theta2.tilde - 
        adjust*cov.mati %*% theta2.var %*% t(cov.mati)/2)/3    
    if(lambda.t2 > 100){
      inf.count = inf.count +1}
    #lambda.t2 <- min(lambda.t2, 0.5) # Set a threshold for now
    Yend <- tk
    dti <- ifelse(deathi == 1, TXi, tk+3)
    dti <- tk+3
    Z2.imp <- 0
    while (Yend <= min(tk+3,dti)) {
      Y2.imp <- rexp(1, lambda.t2)
      Yend <- Yend + Y2.imp
      if(Yend >= min(tk+3,dti)){break}
      event.miss.end <- event.miss.end %>%
        add_row(id=idi, time=Yend, status=2)
      Z2.imp <- Z2.imp + 1
    }
    if(verbose == TRUE){
      cat(glue(
        "Impute {Z2.imp} recurrent events for subject {idi} from group {Ai} ",
        "with rate {round(lambda.t2,4)} ",
        "at ({tk},{tk+3}]"), "\n")
    }
    Z.t3 <- Z.t3 %>%
      rows_update(tibble(id=idi, Z2_t3=Z2.imp), by="id")
    
  }
  
  event.miss.end <- event.miss.end %>%
    arrange(id, time)
  
  ## Impute continuous variable
  Z.df <- Z.t1 %>%
    left_join(
      Z.t2 %>%
        dplyr::select(id, Z1_t2, Z2_t2, Y1_t6, Y2_t6),
      by = "id"
    ) %>%
    left_join(
      Z.t3 %>%
        dplyr::select(id, Z1_t3, Z2_t3, Y1_t9, Y2_t9),
      by = "id"
    ) %>%
    mutate(TX, death)
  
  target.var <- glue("Y1_t{tk+3}")
  
  Z.df.obs <- Z.df %>%
    filter(!is.na(!!sym(target.var))) 
  
  Z.df.mis <- Z.df %>%
    filter(is.na(!!sym(target.var)))
  
  n.obs  <- nrow(Z.df.obs)
  n.mis  <- nrow(Z.df.mis)
  id.obs <- Z.df.obs$id
  id.mis <- Z.df.mis$id
  
  formula3 <- as.formula(
    glue(
      "Y1_t{tk+3}", " ~ ",
      glue_collapse(var.names[[k]]$c,sep = " + ")
    )
  )
  
  fit3 <- lm(formula3, data = Z.df.obs, x=TRUE, y=TRUE)
  theta3.mu <- coef(fit3)
  theta3.var <- vcov(fit3)
  x <- fit3$x
  y <- fit3$y
  nj <- nrow(x)
  pj <- ncol(x)
  L <- chol(solve(t(x) %*% x))
  epsilon.sd <- sqrt(deviance(fit3)/df.residual(fit3))
  
  Zi <- MASS::mvrnorm(1, rep(0,pj), diag(1,nrow=pj))
  gi <- rchisq(1, df=nj-pj-1)
  sigma2i <- epsilon.sd^2*(nj-pj-1)/gi
  sigmai <- sqrt(sigma2i)
  theta3.tilde <- theta3.mu + sigmai*t(L)%*%Zi
  
  cov.mati <- Z.df.mis %>%
    dplyr::select(all_of(var.names[[k]]$c)) %>%
    as.matrix()
  cov.mati <- cbind(1, cov.mati)
  
  Y1_t9.imp <- cov.mati %*% theta3.tilde + rnorm(n.mis)*sigmai
  Z.t3$Y1_t9[Z.t2$id %in% id.mis] <- Y1_t9.imp
  cov.miss$Y1_t9[cov.miss$id %in% id.mis] <- Y1_t9.imp
  
  
  # Y1_t9.imp <- NULL
  # for (i in c(1:n.mis)) {
  #   idi <- id.mis[i]
  #   Zi <- MASS::mvrnorm(1, rep(0,pj), diag(1,nrow=pj))
  #   gi <- rchisq(1, df=nj-pj-1)
  #   sigma2i <- epsilon.sd^2*(nj-pj-1)/gi
  #   sigmai <- sqrt(sigma2i)
  #   theta3.tilde <- theta3.mu + sigmai*t(L)%*%Zi
  #   
  #   cov.mati <- Z.df.mis %>%
  #     filter(id == idi) %>%
  #     dplyr::select(all_of(var.names[[k]]$c)) %>%
  #     as.matrix()
  #   cov.mati <- cbind(1, cov.mati)
  #   Y1_t9.imp <- c(Y1_t9.imp, cov.mati %*% theta3.tilde + rnorm(1)*sigmai)
  # }
  # Z.t3$Y1_t9[Z.t2$id %in% id.mis] <- Y1_t9.imp
  # cov.miss$Y1_t9[cov.miss$id %in% id.mis] <- Y1_t9.imp
  
  # fit3 <- lm(formula3, data = Z.df.obs)
  # theta3.mu <- coef(fit3)
  # theta3.var <- vcov(fit3)
  # 
  # theta3.tilde <-  MASS::mvrnorm(n.mis, theta3.mu, theta3.var)
  # cov.mati <- Z.df.mis %>%
  #   dplyr::select(all_of(var.names[[k]]$c)) %>%
  #   as.matrix()
  # cov.mati <- cbind(1, cov.mati)
  # epsilon.sd <- sqrt(deviance(fit3)/df.residual(fit3))
  # epsilon.Y <- rnorm(n.mis, mean = 0, sd=epsilon.sd)
  # Y1_t9.imp <- apply(cov.mati * theta3.tilde,1,sum) + epsilon.Y
  # Z.t3$Y1_t9[Z.t2$id %in% id.mis] <- Y1_t9.imp
  # cov.miss$Y1_t9[cov.miss$id %in% id.mis] <- Y1_t9.imp
  
  # Impute binary variables
  Z.df <- Z.t1 %>%
    left_join(
      Z.t2 %>%
        dplyr::select(id, Z1_t2, Z2_t2, Y1_t6, Y2_t6),
      by = "id"
    ) %>%
    left_join(
      Z.t3 %>%
        dplyr::select(id, Z1_t3, Z2_t3, Y1_t9, Y2_t9),
      by = "id"
    ) %>%
    mutate(TX, death) 
  
  target.var <- glue("Y2_t{tk+3}")
  
  Z.df.obs <- Z.df %>%
    filter(!is.na(!!sym(target.var))) 
  
  Z.df.mis <- Z.df %>%
    filter(is.na(!!sym(target.var)))
  
  n.obs  <- nrow(Z.df.obs)
  n.mis  <- nrow(Z.df.mis)
  id.obs <- Z.df.obs$id
  id.mis <- Z.df.mis$id
  
  formula4 <- as.formula(
    glue(
      "Y2_t{tk+3}", " ~ ",
      glue_collapse(var.names[[k]]$b,sep = " + ")
    )
  )
  
  fit4 <- glm(formula4, data = Z.df.obs, family = "binomial")
  theta4.mu <- coef(fit4)
  theta4.var <- vcov(fit4)
  
  theta4.tilde <-  MASS::mvrnorm(1, theta4.mu, theta4.var)
  cov.mati <- Z.df.mis %>%
    dplyr::select(all_of(var.names[[k]]$b)) %>%
    as.matrix()
  cov.mati <- cbind(1, cov.mati)
  Y2_t9.imp <- expit(cov.mati %*% theta4.tilde)# + epsilon.Y
  U <- runif(n.mis)
  #Y2_t9.imp <- (Y2_t9.imp>U)*1
  Y2_t9.imp <- (Y2_t9.imp > 0.5)*1
  cov.miss$Y2_t9[cov.miss$id %in% id.mis] <- Y2_t9.imp
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  ### Last interval (9, 12] Impute death time ##
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
  event.miss <- event.miss.end
  k <- 3
  tk <- t[k]
  if(verbose2){cat(glue("Imputation for ({tk}, {tk+3}]"),"\n")}
  
  # Update Z matrix, death time and death status
  Z.mat <- update_Z_mat(event.miss, t)
  Z.lst <- update_Z_lst(cov.miss, Z.mat)
  Z.t1 <- Z.lst$t1
  Z.t2 <- Z.lst$t2
  Z.t3 <- Z.lst$t3
  Z.t4 <- Z.lst$t4
  
  TX <- Z.mat$TX
  death <- Z.mat$death
  
  Z.df <- Z.t1 %>%
    left_join(
      Z.t2 %>%
        dplyr::select(id, Z1_t2, Z2_t2, Y1_t6,Y2_t6),
      by = "id"
    ) %>%
    left_join(
      Z.t3 %>%
        dplyr::select(id, Z1_t3, Z2_t3, Y1_t9,Y2_t9),
      by = "id"
    ) %>%
    left_join(
      Z.t4 %>%
        dplyr::select(id, Z1_t4, Z2_t4, Y1_t12,Y2_t12),
      by = "id"
    ) %>%
    mutate(TX, death) %>%
    filter( !(TX <=tk & death==1))
  
  target.var <- glue("Z1_t{k+1}")
  
  # Z.df.obs <- Z.df %>%
  #   filter(!is.na(!!sym(target.var))) %>%
  #   mutate(time = TX - tk,
  #          status = death)
  
  Z.df.obs <- Z.df %>%
    filter(!is.na(!!sym(target.var)))%>%
    mutate(time = pmin(TX - tk, 3),
           status = ifelse(TX - tk < 3, 1, 0))
  
  Z.df.mis <- Z.df %>%
    filter(is.na(!!sym(target.var)))
  
  n.obs  <- nrow(Z.df.obs)
  id.obs <- Z.df.obs$id
  n.mis  <- nrow(Z.df.mis)
  id.mis <- Z.df.mis$id
  
  formula1 <- as.formula(
    glue(
      "Surv(time, status)", " ~ ",
      glue_collapse(var.names[[k]]$d,sep = " + ")
    )
  )
  
  fit1 <- survreg(
    formula1,
    dist='exponential',
    data = Z.df.obs)
  
  # fit1 <- phreg(
  #   formula1, 
  #   dist = "weibull", shape = 1,
  #   data = Z.df.obs)
  # 
  theta1.mu <- fit1$coefficients
  theta1.var <- fit1$var
  if(random.samp.d){
    theta1.tilde <- MASS::mvrnorm(1, theta1.mu, theta1.var)
  }else{
    theta1.tilde <- theta1.mu 
  }
  event.miss.tmp <- event.miss
  for (i in c(1:n.mis)) {
    idi <- id.mis[i]
    coef.tilde <- theta1.tilde
    cov.mati <- Z.df.mis %>%
      filter(id == idi) %>%
      dplyr::select(all_of(var.names[[k]]$d)) %>%
      as.matrix()
    #cov.mati <- cbind(cov.mati,-1)
    cov.mati <- cbind(-1, -cov.mati)
    ratei <- exp(
      cov.mati %*% coef.tilde - 
        adjust * cov.mati %*% theta1.var %*% t(cov.mati)/2)
    Y1.imp <- rexp(1, rate = ratei)
    Z1.imp <- 1*(Y1.imp <= 3)
    Z.t4 <- Z.t4 %>%
      rows_update(tibble(id=idi, Z1_t4=Z1.imp), by="id")
    if(!is.na(Z1.imp) & Z1.imp == 1){
      event.miss.tmp <- event.miss.tmp %>%
        rows_delete(
          tibble(id=idi,time=tk, status=0),
          by=c("id","time","status")) %>%
        add_row(id=idi, time=Y1.imp+tk, status=1) %>%
        add_row(id=idi, time=tk+3, status=0)
      if(verbose == TRUE){
        cat(glue(
          "Impute {round(Y1.imp+tk,4)} death time with rate {round(ratei,4)} for subject {idi} ",
          "at ({tk},{tk+3}]"), "\n")
      }
    }else{
      event.miss.tmp <- event.miss.tmp %>%
        rows_delete(
          tibble(id=idi,time=tk, status=0),
          by=c("id","time","status")) %>%
        add_row(id=idi, time=tk+3, status=0)
      if(verbose == TRUE){
        cat(glue(
          "Impute {tk+3} censor time with rate {round(ratei,4)} for subject {idi} ",
          "at ({tk},{tk+3}]"), "\n")
      }
      
    }
  }
  event.miss.tmp <- event.miss.tmp %>%
    arrange(id, time)
  
  idzz <- event.miss.tmp$id[event.miss.tmp$time==tk & event.miss.tmp$status==0]
  
  event.miss.tmp <- event.miss.tmp %>%
    rows_delete(tibble(id=idzz,time=tk, status=0),
                by=c("id","time","status")) %>%
    add_row(id=idzz,time=tk+3, status=0) %>%
    arrange(id, time)
  
  ## Impute recurrent time
  ## Update death and censoring time
  Z.mat <- update_Z_mat(event.miss.tmp, t)
  TX <- Z.mat$TX
  death <- Z.mat$death
  
  Z.df <- Z.t1 %>%
    left_join(
      Z.t2 %>%
        dplyr::select(id, Z1_t2, Z2_t2, Y1_t6,Y2_t6),
      by = "id"
    ) %>%
    left_join(
      Z.t3 %>%
        dplyr::select(id, Z1_t3, Z2_t3, Y1_t9,Y2_t9),
      by = "id"
    ) %>% 
    left_join(
      Z.t4 %>%
        dplyr::select(id, Z1_t4, Z2_t4, Y1_t12,Y2_t12),
      by = "id"
    ) %>%
    mutate(TX, death) 
  
  target.var <- glue("Z2_t{k+1}")
  
  Z.df.obs <- Z.df %>%
    filter(!is.na(!!sym(target.var))) 
  
  Z.df.mis <- Z.df %>%
    filter(is.na(!!sym(target.var)))
  
  n.obs <- nrow(Z.df.obs)
  n.mis <- nrow(Z.df.mis)
  id.mis <- Z.df.mis$id
  
  formula2 <- as.formula(
    glue(
      glue("Z2_t{k+1}"), " ~ ",
      glue_collapse(var.names[[k]]$r,sep = " + ")
    )
  )
  
  fit2 <- glm.nb(formula2, data = Z.df.obs)
  theta2.mu <- coef(fit2)
  theta2.var <- vcov(fit2)
  if(random.samp.r){
    theta2.tilde <- MASS::mvrnorm(1, theta2.mu, theta2.var)
  }else{
    theta2.tilde <- theta2.mu
  }
  event.miss.end <- event.miss.tmp
  for (i in c(1:n.mis)) {
    idi <- id.mis[i]
    Ai <- Z.df.mis$A[i]
    TXi <- Z.df.mis$TX[i]
    deathi <- Z.df.mis$death[i]

    cov.mati <- Z.df.mis %>%
      filter(id == idi) %>%
      dplyr::select(all_of(var.names[[k]]$r)) %>%
      as.matrix()
    cov.mati <- cbind(1, cov.mati)
    lambda.t2 <- exp(
      cov.mati %*% theta2.tilde - 
        adjust*cov.mati %*% theta2.var %*% t(cov.mati)/2)/3    
    if(lambda.t2 > 100){inf.count = inf.count +1}
    Yend <- tk
    dti <- ifelse(deathi == 1, TXi, tk+3)
    dti <- tk+3
    Z2.imp <- 0
    while (Yend <= min(tk+3,dti)) {
      Y2.imp <- rexp(1, lambda.t2)
      Yend <- Yend + Y2.imp
      if(Yend >= min(tk+3,dti)){break}
      event.miss.end <- event.miss.end %>%
        add_row(id=idi, time=Yend, status=2)
      Z2.imp <- Z2.imp + 1
    }
    if(verbose == TRUE){
      cat(glue(
        "Impute {Z2.imp} recurrent events for subject {idi} from group{Ai} ",
        "with rate {round(lambda.t2,4)} ",
        "at ({tk},{tk+3}]"), "\n")
    }
    Z.t4 <- Z.t4 %>%
      rows_update(tibble(id=idi, Z2_t4=Z2.imp), by="id")
  }
  event.miss.end <- event.miss.end %>%
    arrange(id, time)
  
  ## Impute continuous variable
  Z.df <- Z.t1 %>%
    left_join(
      Z.t2 %>%
        dplyr::select(id, Z1_t2, Z2_t2, Y1_t6,Y2_t6),
      by = "id"
    ) %>%
    left_join(
      Z.t3 %>%
        dplyr::select(id, Z1_t3, Z2_t3, Y1_t9,Y2_t9),
      by = "id"
    ) %>%
    left_join(
      Z.t4 %>%
        dplyr::select(id, Z1_t4, Z2_t4, Y1_t12,Y2_t12),
      by = "id"
    ) %>%
    mutate(TX, death) 
  
  target.var <- glue("Y1_t{tk+3}")
  
  Z.df.obs <- Z.df %>%
    filter(!is.na(!!sym(target.var))) 
  
  Z.df.mis <- Z.df %>%
    filter(is.na(!!sym(target.var)))
  
  n.obs  <- nrow(Z.df.obs)
  n.mis  <- nrow(Z.df.mis)
  id.obs <- Z.df.obs$id
  id.mis <- Z.df.mis$id
  
  formula3 <- as.formula(
    glue(
      "Y1_t{tk+3}", " ~ ",
      glue_collapse(var.names[[k]]$c,sep = " + ")
    )
  )
  
  fit3 <- lm(formula3, data = Z.df.obs, x=TRUE, y=TRUE)
  theta3.mu <- coef(fit3)
  theta3.var <- vcov(fit3)
  x <- fit3$x
  y <- fit3$y
  nj <- nrow(x)
  pj <- ncol(x)
  L <- chol(solve(t(x) %*% x))
  epsilon.sd <- sqrt(deviance(fit3)/df.residual(fit3))
  
  Zi <- MASS::mvrnorm(1, rep(0,pj), diag(1,nrow=pj))
  gi <- rchisq(1, df=nj-pj-1)
  sigma2i <- epsilon.sd^2*(nj-pj-1)/gi
  sigmai <- sqrt(sigma2i)
  theta3.tilde <- theta3.mu + sigmai*t(L)%*%Zi
  
  cov.mati <- Z.df.mis %>%
    dplyr::select(all_of(var.names[[k]]$c)) %>%
    as.matrix()
  cov.mati <- cbind(1, cov.mati)
  
  Y1_t12.imp <- cov.mati %*% theta3.tilde + rnorm(n.mis)*sigmai
  Z.t4$Y1_t12[Z.t4$id %in% id.mis] <- Y1_t12.imp
  cov.miss$Y1_t12[cov.miss$id %in% id.mis] <- Y1_t12.imp
  
  # Y1_t12.imp <- NULL
  # for (i in c(1:n.mis)) {
  #   idi <- id.mis[i]
  #   Zi <- MASS::mvrnorm(1, rep(0,pj), diag(1,nrow=pj))
  #   gi <- rchisq(1, df=nj-pj-1)
  #   sigma2i <- epsilon.sd^2*(nj-pj-1)/gi
  #   sigmai <- sqrt(sigma2i)
  #   theta3.tilde <- theta3.mu + sigmai*t(L)%*%Zi
  #   
  #   cov.mati <- Z.df.mis %>%
  #     filter(id == idi) %>%
  #     dplyr::select(all_of(var.names[[k]]$c)) %>%
  #     as.matrix()
  #   cov.mati <- cbind(1, cov.mati)
  #   Y1_t12.imp <- c(Y1_t12.imp, cov.mati %*% theta3.tilde + rnorm(1)*sigmai)
  # }
  # Z.t4$Y1_t12[Z.t4$id %in% id.mis] <- Y1_t12.imp
  # cov.miss$Y1_t12[cov.miss$id %in% id.mis] <- Y1_t12.imp
  
  # fit3 <- lm(formula3, data = Z.df.obs)
  # theta3.mu <- coef(fit3)
  # theta3.var <- vcov(fit3)
  # 
  # theta3.tilde <-  MASS::mvrnorm(n.mis, theta3.mu, theta3.var)
  # cov.mati <- Z.df.mis %>%
  #   dplyr::select(all_of(var.names[[k]]$c)) %>%
  #   as.matrix()
  # cov.mati <- cbind(1, cov.mati)
  # epsilon.sd <- sqrt(deviance(fit3)/df.residual(fit3))
  # epsilon.Y <- rnorm(n.mis, mean = 0, sd=epsilon.sd)
  # 
  # Y1_t12.imp <- apply(cov.mati * theta3.tilde,1,sum) + epsilon.Y
  # Z.t4$Y1_t12[Z.t4$id %in% id.mis] <- Y1_t12.imp
  # cov.miss$Y1_t12[cov.miss$id %in% id.mis] <- Y1_t12.imp
  
  # Impute binary variables
  Z.df <- Z.t1 %>%
    left_join(
      Z.t2 %>%
        dplyr::select(id, Z1_t2, Z2_t2, Y1_t6,Y2_t6),
      by = "id"
    ) %>%
    left_join(
      Z.t3 %>%
        dplyr::select(id, Z1_t3, Z2_t3, Y1_t9,Y2_t9),
      by = "id"
    ) %>%
    left_join(
      Z.t4 %>%
        dplyr::select(id, Z1_t4, Z2_t4, Y1_t12,Y2_t12),
      by = "id"
    ) %>%
    mutate(TX, death) 
  
  target.var <- glue("Y2_t{tk+3}")
  
  Z.df.obs <- Z.df %>%
    filter(!is.na(!!sym(target.var))) 
  
  Z.df.mis <- Z.df %>%
    filter(is.na(!!sym(target.var)))
  
  n.obs  <- nrow(Z.df.obs)
  n.mis  <- nrow(Z.df.mis)
  id.obs <- Z.df.obs$id
  id.mis <- Z.df.mis$id
  
  formula4 <- as.formula(
    glue(
      "Y2_t{tk+3}", " ~ ",
      glue_collapse(var.names[[k]]$b,sep = " + ")
    )
  )
  
  fit4 <- glm(formula4, data = Z.df.obs, family = "binomial")
  theta4.mu <- coef(fit4)
  theta4.var <- vcov(fit4)
  
  theta4.tilde <-  MASS::mvrnorm(1, theta4.mu, theta4.var)
  cov.mati <- Z.df.mis %>%
    dplyr::select(all_of(var.names[[k]]$b)) %>%
    as.matrix()
  cov.mati <- cbind(1, cov.mati)
  Y2_t12.imp <- expit(cov.mati %*% theta4.tilde)# + epsilon.Y
  Y2_t12.imp[is.na(Y2_t12.imp)] <- 1
  U <- runif(n.mis)
  #Y2_t12.imp <- (Y2_t12.imp>U)*1
  Y2_t12.imp <- (Y2_t12.imp > 0.5)*1
  cov.miss$Y2_t12[cov.miss$id %in% id.mis] <- Y2_t12.imp
  
  
  event.miss <- event.miss.end
  if(verbose2){cat("end of imputation \n")}
  
  return(list(
    event.imp = event.miss,
    cov.imp = cov.miss,
    inf.count = inf.count  ))
}

mul_imputation <- function(
  event.miss, cov.miss, nimp=50,
  var.names,
  adjust=TRUE, 
  random.samp.d = TRUE, 
  random.samp.r = TRUE,
  verbose=FALSE, 
  verbose2=TRUE){
  
  # id for control and treatment groups
  id.A0 <- cov.data$id[cov.data$A==0]
  id.A1 <- cov.data$id[cov.data$A==1]
  
  event.miss.A0 <- event.miss[event.miss$id %in% id.A0,]
  cov.miss.A0 <- cov.miss[cov.miss$id %in% id.A0, ]
  
  event.miss.A1 <- event.miss[event.miss$id %in% id.A1,]
  cov.miss.A1 <- cov.miss[cov.miss$id %in% id.A1, ]
  
  survest.imp.A0.it <- NULL
  survest.imp.A1.it <- NULL
  survest.dat.A0.it <- NULL
  survest.dat.A1.it <- NULL
  survest.mis.A0.it <- NULL
  survest.mis.A1.it <- NULL
  
  mcfest.imp.A0.it <- NULL
  mcfest.imp.A1.it <- NULL
  mcfest.dat.A0.it <- NULL
  mcfest.dat.A1.it <- NULL
  mcfest.mis.A0.it <- NULL
  mcfest.mis.A1.it <- NULL
  
  Y1.df.it <- NULL
  Y2.df.it <- NULL
  
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
    
    # Summarzing results!!
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
    
    survest.dat.A0.it <- rbind(survest.dat.A0.it, summary(km.dat.A0, times = t.seq)$surv)
    survest.dat.A1.it <- rbind(survest.dat.A1.it, summary(km.dat.A1, times = t.seq)$surv)
    
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
    
    survest.mis.A0.it <- rbind(survest.mis.A0.it, summary(km.mis.A0, times = t.seq)$surv)
    survest.mis.A1.it <- rbind(survest.mis.A1.it, summary(km.mis.A1, times = t.seq)$surv)
    
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
    
    mcfest.imp.A0.it <- rbind(mcfest.imp.A0.it, mcf.stepfun.imp.A0(t.seq))
    mcfest.imp.A1.it <- rbind(mcfest.imp.A1.it, mcf.stepfun.imp.A1(t.seq))
    
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
    
    mcfest.dat.A0.it <- rbind(mcfest.dat.A0.it, mcf.stepfun.dat.A0(t.seq))
    mcfest.dat.A1.it <- rbind(mcfest.dat.A1.it, mcf.stepfun.dat.A1(t.seq))
    
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
    
    mcfest.mis.A0.it <- rbind(mcfest.mis.A0.it, mcf.stepfun.mis.A0(t.seq))
    mcfest.mis.A1.it <- rbind(mcfest.mis.A1.it, mcf.stepfun.mis.A1(t.seq))
    
    Y1.df.i <- cov.miss %>%
      mutate(type = "No imputation") %>%
      bind_rows(cov.data %>% mutate(type="No censoring")) %>%
      bind_rows(cov.imp %>% mutate(type="Imputation")) %>%
      group_by(type, A) %>%
      summarise(
        across(starts_with("Y1"), ~mean(.x, na.rm=TRUE)),
        .groups = 'drop') %>%
      mutate(iteration = i) %>%
      pivot_longer(
        starts_with("Y"), names_prefix = "Y1_t",
        names_to = "time", values_to = "Y1")
    
    Y1.df.it <- rbind.data.frame(Y1.df.it, Y1.df.i)
    
    Y2.df.i <- cov.miss %>%
      mutate(type = "No imputation") %>%
      bind_rows(cov.data %>% mutate(type="No censoring")) %>%
      bind_rows(cov.imp %>% mutate(type="Imputation")) %>%
      group_by(type, A) %>%
      summarise(
        across(starts_with("Y2"), ~mean(.x, na.rm=TRUE)),
        .groups = 'drop') %>%
      mutate(iteration = i) %>%
      pivot_longer(
        starts_with("Y"), names_prefix = "Y2_t",
        names_to = "time", values_to = "Y2")
    
    Y2.df.it <- rbind.data.frame(Y2.df.it, Y2.df.i)
  }
  
  return(list(
    survest.imp.A0.it = survest.imp.A0.it,
    survest.imp.A1.it = survest.imp.A1.it,
    survest.dat.A0.it = survest.dat.A0.it,
    survest.dat.A1.it = survest.dat.A1.it,
    survest.mis.A0.it = survest.mis.A0.it,
    survest.mis.A1.it = survest.mis.A1.it,
    
    mcfest.imp.A0.it = mcfest.imp.A0.it,
    mcfest.imp.A1.it = mcfest.imp.A1.it,
    mcfest.dat.A0.it = mcfest.dat.A0.it,
    mcfest.dat.A1.it = mcfest.dat.A1.it,
    mcfest.mis.A0.it = mcfest.mis.A0.it,
    mcfest.mis.A1.it = mcfest.mis.A1.it,
    
    Y1.df.it = Y1.df.it,
    Y2.df.it = Y2.df.it
  ))
  
}



