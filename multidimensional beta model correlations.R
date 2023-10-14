nimble_functions <- function(){
  
  require("nimble")
  
  # multivariate Gaussian copula (density)
  dcopula_normal <- nimbleFunction(
    run = function(target=double(1),par = character(0, default="par"),u = double(1), SIGMA = double(2), prec_param = integer(0, default = 0),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      if (par=="u"){
        u = target
      }else if (par=="SIGMA"){
        SIGMA[lower.tri(SIGMA,diag=T)] = target
        SIGMA[upper.tri(SIGMA)] = SIGMA[lower.tri(SIGMA)]
      }
      # u is a random variable within [0,1]
      x = qnorm(u)
      term1 = sum(dnorm(x, log=TRUE))
      cholesky = chol(SIGMA)
      term2 = dmnorm_chol(x = x, mean = rep(0, dim(cholesky)[1]), cholesky = cholesky, prec_param = prec_param, log = TRUE)
      logp = term2 - term1
      if(log) return(logp)
      else return(exp(logp))
    })
  # Gaussian copula (cumulative distribution)
  pcopula_normal <- nimbleFunction(
    run = function(u = double(1), mu = double(1), SIGMA = double(2)) {
      returnType(double(0))
      # qnorm function: Transforming (0,1) to normal scale (-Inf Inf)
      x <- qnorm(u)
      return(mvtnorm::pmvnorm(upper=x, mean = mu, sigma = SIGMA)[1])
  })
  
  
  # test
  # u = c(0.4,0.6) # u is a random variable within [0,1]
  u = c(0.3667176, 0.2140250)
  MU = rep(0,length(u))
  corr = 0.5
  SIGMA = matrix(corr,length(u),length(u))
  diag(SIGMA) = 1
  dcopula_normal(target=u,par="u",u=u, SIGMA = SIGMA, log=FALSE)
  pcopula_normal(u, MU, SIGMA)
  
  # NOTE:
  # upper-triangular Cholesky factor of either the precision matrix (when prec_param is TRUE) or covariance matrix (otherwise).

  
  # density: multivariate beta copula density
  dmbeta <- nimbleFunction(
    run = function(x = double(1), a = double(1), b = double(1), cholesky = double(2), prec_param = integer(0, default = 0),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      # x is the random variable from beta density
      u = pbeta(x, a, b)
      qq = qnorm(u)
      qq[qq == Inf] = 10 # in case Inf
      term1 = sum(dnorm(qq, log=TRUE))
      term2 = dmnorm_chol(x = qq, mean = rep(0, dim(cholesky)[1]), cholesky = cholesky, prec_param = prec_param, log = TRUE)
      logp = sum(dbeta(x, a, b, log=TRUE)) + term2 - term1
      # if (logp < -500){logp = -500}
      if(log) return(logp)
      else return(exp(logp))
    })
  # cumulative distribution: multivariate beta distribution
  pmbeta <- function(x = double(1), a = double(1), b = double(1), SIGMA = double(2)){
      # qbeta function: Transforming (0,1) to normal scale (-Inf Inf)
      u <- pbeta(x, a, b)
      qq = qnorm(u)
      return(mvtnorm::pmvnorm(upper=qq, mean = rep(0, dim(SIGMA)[1]), sigma = SIGMA)[1])
  }
  
  # test
  x = c(0.4,0.6) # b is the random variable from beta density
  MU = rep(0,length(x))
  corr = 0.5
  SIGMA = matrix(corr,length(x),length(x))
  diag(SIGMA) = 1
  a = c(2, 2)
  b = c(5, 2)
  dmbeta(x=x, a=a, b=b, cholesky = chol(SIGMA), prec_param=0, log=FALSE)
  pmbeta(x, a, b, SIGMA)
  
  assign('dmbeta', dmbeta, envir = .GlobalEnv)
  assign('pmbeta', pmbeta, envir = .GlobalEnv)
 
  
  
  # random number: multivariate beta copula distribution
  rmbeta <- nimbleFunction(
    run = function(n = integer(0, default = 1), a = double(1), b = double(1), cholesky = double(2), prec_param = integer(0, default = 0)) {
      returnType(double(1))
      z <- rmnorm_chol(n = n, mean=rep(0, dim(cholesky)[1]), cholesky = cholesky, prec_param = prec_param)
      COV <- t(cholesky)%*%cholesky
      return(qbeta(pnorm(z, 0, sqrt(diag(COV))), a, b))
    })
  # test
  x = c(0.4,0.6)
  MU = rep(0,length(x))
  corr = 0
  SIGMA = matrix(corr,length(x),length(x))
  diag(SIGMA) = 1
  a = c(2, 2)
  b = c(5, 2)
  rmbeta(n=1, a=a, b=b, cholesky = chol(SIGMA))
  
  assign('rmbeta', rmbeta, envir = .GlobalEnv)
 
  
 
  
  
  ##################### in nimble
  # density: multivariate logit-normal distribution
  dlogitnormal <- nimbleFunction(
    run = function(x = double(1), mu = double(1), cholesky = double(2), prec_param = integer(0, default = 0),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      logx = log(x/(1-x))
      # logx = qlogis(x) # same with log(x/(1-x))
      logp = - sum(log(x*(1-x))) + dmnorm_chol(x = logx, mean = mu, cholesky = cholesky, prec_param = prec_param, log = TRUE)
      if(log) return(logp)
      else return(exp(logp))
    })
  # test
  x = c(0.4,0.6)
  MU = rep(0,length(x))
  MU[1] <- -1
  corr = -0.9
  SIGMA = matrix(corr,length(x),length(x))
  SIGMA[1,1] <- 1.5^2
  SIGMA[2,2] <- 1^2
  dlogitnormal(x, mu = MU, cholesky = chol(SIGMA), log=FALSE)
  
  assign('dlogitnormal', dlogitnormal, envir = .GlobalEnv)
  
  
  # (VERSION 1) density: multivariate logit-normal copula density
  # Use Gaussian copula to link many univariate logit-normal distribution with assigned correlations!
  dlogitnormal_copula1 <- nimbleFunction(
    run = function(x = double(1), mu = double(1), SIGMA = double(2), prec_param = integer(0, default = 0),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      # x is the random variable from beta density
      logit_x <- log(x/(1-x))
      var = diag(SIGMA)
      v <- (logit_x-mu)/sqrt(2*var)
      erf_v <- pchisq(2 * v^2, 1) * sign(v) # same with pracma::erf(v)
      u = 0.5*(1 + erf_v)
      u_zero = which(u==0)
      u_one = which(u==1)
      if (length(u_zero)!=0) u[u_zero] = 1e-8
      if (length(u_one)!=0) u[u_one] = 1 - 1e-8
      qq = qnorm(u)
      qq[qq == Inf] = 10 # in case Inf
      term1 = sum(dnorm(qq, log=TRUE))
      SIGMA <- cov2cor(SIGMA) # IMPORTANT: If var is set, we must transform the covariance matrix to correlation matrix
      cholesky = chol(SIGMA)
      term2 = dmnorm_chol(x = qq, mean = rep(0, dim(cholesky)[1]), cholesky = cholesky, prec_param = prec_param, log = TRUE)
      
      logp_logit_normal = - log(x*(1-x))  + dnorm(x = logit_x, mean = mu, sd = sqrt(var), log = TRUE)
      
      logp = sum(logp_logit_normal) + term2 - term1
      # if (logp < -500){logp = -500}
      if(log) return(logp)
      else return(exp(logp))
    })
  # test
  dlogitnormal_copula1(x=x, mu = MU, SIGMA = SIGMA, prec_param=0, log=FALSE)
  
  assign('dlogitnormal_copula1', dlogitnormal_copula1, envir = .GlobalEnv)
  
  
  
  # (VERSION 2) density: multivariate logit-normal copula density
  # Use Gaussian copula to link many univariate logit-normal distribution with assigned correlations!
  dlogitnormal_copula2 <- nimbleFunction(
    run = function(x = double(1), mu = double(1), SIGMA = double(2), prec_param = integer(0, default = 0),
                   log = integer(0, default = 0)) {
      returnType(double(0))
      # x is the random variable from beta density
      logit_x <- log(x/(1-x))
      var = 1 # IMPORTANT: if var is set to 1, the covariance matrix (SIGMA) can be used directly.
      v <- (logit_x-mu)/sqrt(2*var)
      erf_v <- pchisq(2 * v^2, 1) * sign(v) # same with pracma::erf(v)
      u = 0.5*(1 + erf_v)
      u_zero = which(u==0)
      u_one = which(u==1)
      if (length(u_zero)!=0) u[u_zero] = 1e-8
      if (length(u_one)!=0) u[u_one] = 1 - 1e-8
      qq = qnorm(u)
      qq[qq == Inf] = 10 # in case Inf
      term1 = sum(dnorm(qq, log=TRUE))
      cholesky = chol(SIGMA)
      term2 = dmnorm_chol(x = qq, mean = rep(0, dim(cholesky)[1]), cholesky = cholesky, prec_param = prec_param, log = TRUE)
      
      logp_logit_normal = - log(x*(1-x))  + dnorm(x = logit_x, mean = mu, sd = sqrt(var), log = TRUE)
      
      logp = sum(logp_logit_normal) + term2 - term1
      # if (is.nan(exp(logp))) browser()
      if(log) return(logp)
      else return(exp(logp))
    })
  # test
  dlogitnormal_copula2(x=x, mu = MU, SIGMA = SIGMA, prec_param=0, log=FALSE)
  
  assign('dlogitnormal_copula2', dlogitnormal_copula2, envir = .GlobalEnv)
  
  
  # Univariate cumulaitve distribution of logitnormal density
  plogitnormal_univariate <- function(x,mu,var){
    logit_x <- log(x/(1-x))
    v <- (logit_x-mu)/sqrt(2*var)
    erf_v <- pchisq(2 * v^2, 1) * sign(v) # same with pracma::erf(v)
    u = 0.5*(1 + erf_v)
    return(u)
  }
  assign('plogitnormal_univariate', plogitnormal_univariate, envir = .GlobalEnv)
  
  
  # multivariate cumulative logit-normal distribution
  plogitnormal_multivariate <- function(x, mu, sigma){
    # logit function: Transforming (0,1) to normal scale (-Inf Inf)
    qq <- qlogis(x) # same with log(x/(1-x))
    return(mvtnorm::pmvnorm(upper=qq, mean = mu, sigma = sigma)[1])
  }
  plogitnormal_multivariate(x, MU, SIGMA)
  
  assign('plogitnormal_multivariate', plogitnormal_multivariate, envir = .GlobalEnv)
  
  # # (numerical integral) IMPORTANT: The first argument of function must not be included when using adaptIntegrate!!
  # # Check if the integral of logit-normal density is equal to the cumulative logit-normal distribution
  # # Result is equal to plogitnormal_multivariate(x, MU, SIGMA) above
  # library("cubature")
  # adaptIntegrate(dlogitnormal, lowerLimit = c(0, 0), upperLimit = x , tol=1e-7, mu = MU, cholesky = chol(SIGMA), log=FALSE)
  # adaptIntegrate(dlogitnormal_copula1, lowerLimit = c(0, 0), upperLimit = x , tol=1e-7, mu = MU, SIGMA = SIGMA, log=FALSE)
  # cubintegrate(dlogitnormal_copula2, lower = c(0, 0), upper = x , relTol=1e-7, method = c("hcubature"), mu = MU, SIGMA = SIGMA, log=FALSE)
  
  
  # multivariate quantile logit-normal distribution
  qlogitnormal_multivariate <- function(x, mu, sigma){
    QQ = mvtnorm::qmvnorm(p=x, mean = mu, sigma = sigma)$quantile
    return(plogis(QQ)) # formula is 1/(1 + exp(-u))
  }
  assign('qlogitnormal_multivariate', qlogitnormal_multivariate, envir = .GlobalEnv)
  
  
  
  
  
  
  # random number: multivariate logit-normal distribution
  rlogitnormal <- nimbleFunction(
    run = function(n = integer(0), mu = double(1), cholesky = double(2), prec_param = integer(0, default = 0)) {
      returnType(double(1))
      return(plogis(rmnorm_chol(n = n, mean=mu, cholesky = cholesky, prec_param = prec_param)))
      # plogis is 1/(1 + exp(-u))
    })
  # test
  rlogitnormal(n=1, MU, cholesky = chol(SIGMA))
  
  assign('rlogitnormal', rlogitnormal, envir = .GlobalEnv)
  
  
  
  
  uppertri_mult_diag <- nimbleFunction(
    run = function(mat = double(2), vec = double(1)) {
      returnType(double(2))
      p <- length(vec)
      out <- matrix(nrow = p, ncol = p, init = FALSE)
      for(i in 1:p)
        out[ , i] <- mat[ , i] * vec[i]
      return(out)
    })
  
  assign('uppertri_mult_diag', uppertri_mult_diag, envir = .GlobalEnv)
  
  
}

Conditional_Predictive_Ordinate_fit <- function(prob){
  
  L = dim(prob)[1]
  w <- 1/prob
  w_bar <- apply(w, c(2,3), mean)
  w_wave <- pmin(w, sqrt(L)*w_bar) # L is the chain length
  # cpo_2014 is an optimized version of cpo_origin; # see paper: Bayesian beta regression for bounded responses with unknown supports
  cpo_2014 = apply(prob*w_wave, c(2,3), sum) / apply(w_wave, c(2,3), sum) 
  cpo_origin = 1/w_bar
  
  # Larger value of logarithm of the Pseudo-marginal likelihood (LPML) suggests a better fit:
  LPML_2014 = sum(log(cpo_2014)) 
  LPML_origin = sum(log(cpo_origin)) # contain 0s, so log(0) = -Inf
  
  return(list(cpo_2014=cpo_2014
              ,cpo_origin=cpo_origin
              ,LPML_2014=LPML_2014
              ,LPML_origin=LPML_origin))
}

Cox_Snell_fit <- function(t
                          ,r_Cox_Snell
                          ,eap_r_Cox_Snell){
  L = dim(r_Cox_Snell)[1]
  people = dim(r_Cox_Snell)[2]
  total_items = dim(r_Cox_Snell)[3]

  t_L = length(t)
  lambda_item <- array(NA,c(t_L,L,total_items))
  lambda_person <- array(NA,c(t_L,L,people))
  
  lambda_item_eap <- array(NA,c(t_L,total_items))
  lambda_person_eap <- array(NA,c(t_L,people))
  for (k in 1:t_L){
    comp <- r_Cox_Snell <= t[k]
    lambda_item[k,,] = -log(1 - apply(comp,c(1,3),mean))
    lambda_person[k,,] = -log(1 - apply(comp,c(1,2),mean))
    
    # for eap
    comp_eap <- eap_r_Cox_Snell <= t[k]
    lambda_item_eap[k,] = -log(1 - apply(comp_eap,c(2),mean))
    lambda_person_eap[k,] = -log(1 - apply(comp_eap,c(1),mean))
  }
  
  lambda_item[is.infinite(lambda_item)] = NA
  lambda_person[is.infinite(lambda_person)] = NA
  lambda_item_eap[is.infinite(lambda_item_eap)] = NA
  lambda_person_eap[is.infinite(lambda_person_eap)] = NA
  
  return(list(t=t
              ,lambda_item=lambda_item
              ,lambda_person=lambda_person
              ,lambda_item_eap=lambda_item_eap
              ,lambda_person_eap=lambda_person_eap))
}

# calculate M and N
beta_M_N <- function(obs, model, item_index, samples_matrix, keep,people,itemnum,testlet,total_items){
  
  # load customized distributions
  nimble_functions()
  
  D = itemnum
  
  # dim(rep) <- c(dim(rep)[1],people,total_items)
  
  L <- dim(samples_matrix)[1]
  name <- colnames(samples_matrix)
  
  ind <- grepl("alpha", name)
  alpha <- samples_matrix[,ind]
  dim(alpha) <- c(L,total_items,itemnum)
  
  ind <- grepl("theta", name)
  theta <- samples_matrix[,ind]
  dim(theta) <- c(L,people,D)
  
  ind <- grepl("delta", name)
  delta <- samples_matrix[,ind]
  
  ind <- grepl("tau", name)
  tau <- samples_matrix[,ind]
  
  sigma_item_cor <- sigma_item_cor_4D <- matrix(0) # for model == 1
  if (model==0 | model == 2){
    ind <- grepl("sigma_item_cor", name)
    sigma_item_cor <- samples_matrix[,ind]
    sigma_item_cor_4D <- sigma_item_cor
    dim(sigma_item_cor_4D) <- c(L,itemnum,itemnum,testlet)
    
    # obtain EAP (used for Cox-Snell plots)
    eap_sigma_item_cor_4D = apply(sigma_item_cor_4D,c(2,3,4),mean)
  }else{
    eap_sigma_item_cor_4D = NULL
  }
  
  # obtain EAP (used for Cox-Snell plots)
  eap_alpha = apply(alpha,c(2,3),mean)
  eap_theta = apply(theta,c(2,3),mean)
  eap_delta = apply(delta,2,mean)
  eap_tau = apply(tau,2,mean)
  
  M_eap <- N_eap <- eap_r_Cox_Snell <- array(NA,c(people,total_items))
  cov_obs <- cov_rep <- NULL
  
  
  
  
  

  
  # rcpp
  print("Compiling Rcpp code...")
  Rcpp::sourceCpp("sim_beta.cpp")
  
  print("Computing outfit and infit...")
  times <- system.time(out <- sim_beta(obs, model, L, people, itemnum, testlet, total_items, alpha, theta, delta, tau, sigma_item_cor))
  
  
  # # Conditional Predictive Ordinate (CPO); see page 344, 375
  # # Defined as Pr(y_i| y_(i)), where y_(i) means i_th data point is omitted. 
  # Conditional_Predictive_Ordinate_fit(out$prob)
  
  
  # for Cox-Snell plots
  if (model == 0){
    
    for(n in 1:people){
      cnt = 1
      for (k in 1:testlet){
        var <- diag(eap_sigma_item_cor_4D[,,k])
        for (i in 1:itemnum){
          x = obs[n,cnt]
          
          mu <- as.vector(eap_delta[cnt] + eap_alpha[cnt,]%*%eap_theta[n,])
          u = plogitnormal_univariate(x,mu,var[i])
 
          # for Cox-Snell plots
          eap_r_Cox_Snell[n,cnt] <- -log(1 - u)
          
          cnt = cnt + 1
        }
      }
    }
    
  }else if (model==1 | model == 2){
    for(n in 1:people){
      for (k in 1:total_items){
        M_eap[n,k] <- exp(((sum(eap_alpha[k,1:D]*eap_theta[n,1:D])-eap_delta[k])+eap_tau[k])/2)
        N_eap[n,k] <- exp((-(sum(eap_alpha[k,1:D]*eap_theta[n,1:D])-eap_delta[k])+eap_tau[k])/2)
        eap_r_Cox_Snell[n,k] <- -log(1 - pbeta(obs[n,k], M_eap[n,k], N_eap[n,k]))
      }
    }
  }
  

  
  return(list(obs = obs
              , r_rep = out$r_rep
              , nu_obs2=out$nu_obs2
              , nu_rep2=out$nu_rep2
              , var=out$var
              , times=times
              , cov_obs=out$cov_obs
              , cov_rep=out$cov_rep
              , sgddm_obs=out$sgddm_obs
              , sgddm_rep=out$sgddm_rep
         ,eap_sigma_item_cor_4D=eap_sigma_item_cor_4D
         ,eap_alpha=eap_alpha
         ,eap_theta=eap_theta
         ,eap_delta=eap_delta
         ,eap_tau=eap_tau
         
         , prob = out$prob
         , r_Cox_Snell = out$r_Cox_Snell
         , eap_r_Cox_Snell = eap_r_Cox_Snell
         , pred_lower_yi = out$pred_lower_yi
         
         # ,st=st
         ))
  

}

mcmc_main <- function(response=NULL
                      ,rand.seed = 1
                      ,testlet = NULL
                      ,itemnum = 6
                      ,people = 1000
                      ,simulation = FALSE
                      ,return_sim=FALSE
                      ,burnin = 1
                      ,keep = 2
                      ,nchains = 1
                      ,model = 1
                      ,est_theta = FALSE
                      ,sim_data = FALSE
                      ,use_real_data_pars = -1
                      ,use_real_data_pars_prior = FALSE
                      ,get_item_fit=FALSE
                      ,waic_nimble=FALSE
                      ,sigma_item_cor_values = FALSE
                      ,mcmc_length_add = 0){
  
  # load customized distributions
  nimble_functions()
  assign('beta_M_N', beta_M_N, envir = .GlobalEnv)
  
  if (itemnum == 1){
    stop("At least two items in a testlet for current program, because dmbeta and others need to be modified...")
  }
  
  
  if (!is.null(response)){
    people = nrow(response)
    col = ncol(response)
    if (!is.null(testlet) & testlet%%1==0){
      itemnum = col/testlet
    }else{
      stop("Please input a integer 'testlet'. The item number in each testlet must be the same.")
    }
    if (itemnum%%1==0){
      
    }else{
      stop("The item number in each testlet must be the same.")
    }
    D <- itemnum
    total_items = testlet * itemnum
    y = qlogis(response) # same with log(x/(1-x))
    
  }
  
  # simulation == FALSE means analyzing real data
  if (is.null(response) & simulation == FALSE){
    print("Read your data from disk: File name must be 'y01_real.rds'.")
    
    # real data: career interest
    response = readRDS(file = "y01_real.rds")
    

    testlet = 5

    people = nrow(response)
    col = ncol(response)
    if (!is.null(testlet) & testlet%%1==0){
      itemnum = col/testlet # equal to dimension of theta
    }else{
      stop("Please input a integer 'testlet'. The item number in each testlet must be the same.")
    }
    if (itemnum%%1==0){
      
    }else{
      stop("The item number in each testlet must be the same.")
    }

    D <- itemnum
    people <- nrow(response)
    
    total_items = testlet * itemnum
    
    # logit function: Transforming (0,1) to normal scale (-Inf Inf)
    y = qlogis(response) # same with log(x/(1-x))
  }
  

  # use_real_data_pars != -1 means to use item estimates from real data analysis
  if (use_real_data_pars != -1){
    itemnum <- 6 # equal to dimension of theta
    D <- itemnum
    testlet <- 5
    total_items = testlet * itemnum
    
    if (use_real_data_pars == TRUE){
      stop("Set use_real_data_pars = -1 if use_real_data_pars is TRUE.")
    }
  }
  
  set.seed(1)

  # itemnum <- 6 # equal to dimension of theta
  D <- itemnum
  
  # people <- 1000
  
  # testlet <- 5
  
  total_items = testlet * itemnum
  
  library("mvtnorm")
  library("randcorr")
  
  
  # LKJ distribution
  eta = 1
  
  
  
  # If item correlation is considered:
  
  # model = 1 # multidimensional beta model (no correlation)
  # model = 2 # beta copula with item correlations
  # model = 3 # beta copula with no item correlations
  # model = 4 # correlated traits–correlated uniqueness model (CTCUM; Marsh, 1989)

  if (model == 0){
    model_is = "logit-normal model"
  }else if (model == 1){
    model_is = "multidimensional beta model (no correlation)"
  }else if (model == 2){
    model_is = "beta copula with item correlations"
  }else if (model == 3){
    model_is = "beta copula with no item correlations"
  }else if (model == 4){
    model_is = "CTCUM"
  }
  print(model_is)
  

  
  
  


  
  # item error correlations for each testlet
  L22 = array(NA,c(itemnum,itemnum,testlet))
  sigma_item_cor = sigma_item_chol = vector("list",testlet)
  for (k in 1:testlet){
    sigma_item_cor[[k]] = matrix( 0.0 ,itemnum,itemnum)
    if (D==2){
      sigma_item_cor[[k]][] = 0.5
      diag(sigma_item_cor[[k]]) = 1
    }else{
      sigma_item_cor[[k]] = randcorr(D)
    }
    
    if (model == 1 | model == 3){
      sigma_item_cor[[k]] = diag(D)
    }
    
    if ( isTRUE(sigma_item_cor_values) ){
      sigma_item_cor[[k]] = diag(D)
    }
    
    # ONLY work for two items in each testlet
    if ( is.numeric(sigma_item_cor_values) ){
      if (model == 0 |model == 2){
        sigma_item_cor[[k]][upper.tri(sigma_item_cor[[k]])] = sigma_item_cor_values
        sigma_item_cor[[k]][lower.tri(sigma_item_cor[[k]])] = sigma_item_cor_values
      }
    }
    
    sigma_item_chol[[k]] = chol(sigma_item_cor[[k]])
    L22[,,k] = diag(D) # starting value for cholesky matrix
  }
  
  
  
  
  # prior for item parameters (from real data) for simulation
  prior_alpha = prior_delta = prior_tau = prior_sigma_item_cor = list()
  
  prior_alpha$mu = 0
  prior_alpha$sd = 0.3
  prior_delta$mu = 0
  prior_delta$sd = 0.3
  prior_tau$mu = 0
  prior_tau$sd = 0.3
  
  if (simulation == TRUE & use_real_data_pars_prior == TRUE & model!=4){
    
    if ( (model==2 | model==1 | model==3) & file.exists("use_real_data_pars_prior_BCM.rds")){
      real_est = readRDS("use_real_data_pars_prior_BCM.rds")
    } 
    if (model==0 & file.exists("use_real_data_pars_prior_LNM.rds")){
      real_est = readRDS("use_real_data_pars_prior_LNM.rds")
    } 
    
    temp = numeric()
    for (k in 1:testlet){
      temp2 = real_est$alpha[[k]]
      temp = c(temp, temp2[temp2!=0])
    }
    
    library("MASS")
    f1 <- fitdistr(temp,"lognormal")
    prior_alpha$mu = f1$estimate[1]
    prior_alpha$sd = f1$estimate[2]
    
    temp = do.call(c,real_est$delta)
    f1 <- fitdistr(temp,"normal")
    prior_delta$mu = f1$estimate[1]
    prior_delta$sd = f1$estimate[2]
    
    if (model==2){
      temp = do.call(c,real_est$tau)
      f1 <- fitdistr(temp,"normal")
      prior_tau$mu = f1$estimate[1]
      prior_tau$sd = f1$estimate[2]
    }

    
  }else{
    print("no file: use_real_data_pars_prior.rds")
  }
  
  
  
  
  if ( (use_real_data_pars == 0 | model == 0 | model == 4) & use_real_data_pars <= 0 ){
    # item parameters
    alpha = vector("list",testlet)
    design = vector("list",testlet)
    DESIGN = numeric()
    delta = vector("list",testlet)
    # sigma_item_cor = vector("list",testlet)
    error = vector("list",testlet)
    for (k in 1:testlet){
      design[[k]] = diag(D)
      DESIGN = rbind(DESIGN,design[[k]])
      alpha[[k]] = design[[k]] * rlnorm(n=itemnum,prior_alpha$mu,prior_alpha$sd)
      delta[[k]] = rnorm(n=itemnum,prior_delta$mu,prior_delta$sd)
      # sigma_item_cor[[k]] = randcorr(D)
      error[[k]] <- rmvnorm(n=people, mean=rep(0,D), sigma=sigma_item_cor[[k]])
    }
    

    
    # # example:
    # x = c(0.4,0.6)
    # 
    # # logit function: Transforming (0,1) to normal scale (-Inf Inf)
    # u = qlogis(x) # same with log(x/(1-x))
    # 
    # # inverse logit function: Transforming (-Inf,Inf) to original scale (0,1)
    # plogis(u) # equal to x; formula is 1/(1 + exp(-u))
    
    
  }else if (use_real_data_pars == 1 | use_real_data_pars == 2 | use_real_data_pars == 3 
            | model == 1 | model == 2 | model == 3){
    # item parameters
    alpha = vector("list",testlet)
    design = vector("list",testlet)
    DESIGN = numeric()
    delta = vector("list",testlet)
    tau = vector("list",testlet)
    for (k in 1:testlet){
      design[[k]] = diag(D)
      DESIGN = rbind(DESIGN,design[[k]])
      alpha[[k]] = design[[k]] * matrix(rlnorm(n=itemnum*D,prior_alpha$mu,prior_alpha$sd),itemnum,D)
      delta[[k]] = rep(NA,itemnum)
      tau[[k]] = rep(NA,itemnum)
      for (i in 1:itemnum){
        delta[[k]][i] <- rnorm(n=1,prior_delta$mu,prior_delta$sd)
        tau[[k]][i] <- rnorm(n=1,prior_tau$mu,prior_tau$sd)
      }
    }
    

    
    # alpha starting values
    alpha_ini = do.call(rbind,alpha)
    alpha_ini[alpha_ini == 0] = NA
    alpha_ini[alpha_ini != 0] = 1
  }
  

  # same code like above (used for use_real_data_pars==0 and model!=0)
  # for beta models
  design = vector("list",testlet)
  DESIGN = numeric()
  for (k in 1:testlet){
    design[[k]] = diag(D)
    DESIGN = rbind(DESIGN,design[[k]])
  }
  # non-zero slope index
  index_alpha = which(DESIGN==1, arr.ind=TRUE)
  index_zero = which(DESIGN==0, arr.ind=TRUE)
  
  # zero slope index
  index_alpha_length = dim(index_alpha)[1]
  index_zero_length = dim(index_zero)[1]
  
  
  
  # for logit-normal model
  Ustar = array(NA,c(itemnum,itemnum,testlet))
  for (k in 1:testlet){
    Ustar[,,k] = diag(D) # starting value for cholesky matrix
  }
  
  
  
  
  # use the item estimates from real data analysis
  if (use_real_data_pars != -1){
    if (use_real_data_pars == 0){
      est_real = readRDS("simulation_FALSE_model_0_testlet_5_itemnum_6_people_550.rds")
    }else if (use_real_data_pars == 1){
      est_real = readRDS("simulation_FALSE_model_1_testlet_5_itemnum_6_people_550.rds")
    }else if (use_real_data_pars == 2){
      temp = readRDS("simulation_FALSE_model_2_testlet_5_itemnum_6_people_550_item_estimates.rds")
      est_real = list()
      est_real$est_item$statistics = temp
    }else if (use_real_data_pars == 3){
      est_real = readRDS("simulation_FALSE_model_3_testlet_5_itemnum_6_people_550.rds")
    }
    
    est_real_all = round( est_real$est_item$statistics[,"Mean"] , digits = 2)
    est_real_sd = round( est_real$est_item$statistics[,"SD"] , digits = 2)
    
    # est
    ind <- grepl('alpha', names(est_real_all))
    est_real_alpha = est_real_all[ind]
    alpha_se = alpha
    est_real_alpha_sd = est_real_sd[ind]
    
    if (use_real_data_pars == 0){
      cnt = 1:itemnum
      for (k in 1:testlet){
        diag(alpha[[k]]) = est_real_alpha[cnt]
        diag(alpha_se[[k]]) = est_real_alpha_sd[cnt]
        cnt = cnt + itemnum
      }
    }
    
    if (use_real_data_pars == 2){
      cnt = 1:itemnum
      temp = matrix(est_real_alpha, total_items, itemnum)
      for (k in 1:testlet){
        alpha[[k]] = temp[cnt,]
        cnt = cnt + itemnum
      }
    }

    
    # est
    ind <- grepl('delta', names(est_real_all))
    est_real_delta = est_real_all[ind]
    delta = relist( est_real_delta , delta)
    # SE
    delta_se = delta
    est_real_delta = est_real_sd[ind]
    delta_se = relist( est_real_delta , delta_se)
    
    

    if (use_real_data_pars == 1 | use_real_data_pars == 2 | use_real_data_pars == 3){
      # est
      ind <- grepl('tau', names(est_real_all))
      est_real_tau = est_real_all[ind]
      tau = relist( est_real_tau , tau)
      # SE
      tau_se = tau
      est_real_tau = est_real_sd[ind]
      tau_se = relist( est_real_tau , tau_se)
    }else{
      tau = NULL
      tau_se = NULL
    }

    
    # est
    ind <- grepl('sigma_item_cor', names(est_real_all))
    est_real_sigma_item_cor = est_real_all[ind]
    est_real_sigma_item_cor = matrix(est_real_sigma_item_cor, itemnum^2, testlet)
    sigma_item_cor = relist( est_real_sigma_item_cor , sigma_item_cor)
    # SE
    sigma_item_cor_se = sigma_item_cor
    est_real_sigma_item_cor = est_real_sd[ind]
    est_real_sigma_item_cor = matrix(est_real_sigma_item_cor, itemnum^2, testlet)
    sigma_item_cor_se = relist( est_real_sigma_item_cor , sigma_item_cor_se)
    
    
    if (use_real_data_pars == 1 | use_real_data_pars == 3){
      for (k in 1:length(sigma_item_cor)){
        sigma_item_cor[[k]] = diag(itemnum)
      }
    }
    
    
    if (use_real_data_pars == 0){
      error = vector("list",testlet)
      for (k in 1:testlet){
        error[[k]] <- rmvnorm(n=people, mean=rep(0,D), sigma=sigma_item_cor[[k]])
      }
    }
    
    
    # # save estimates for simulation setting
    # saveRDS(list(alpha=alpha,delta=delta,tau=tau,sigma_item_cor=sigma_item_cor),file=paste("use_real_data_pars_prior",".rds", sep = ""))
    
    
  }


  
  
  
  

  
  
  
  
  # latent traits

  mu = rep(0,D)
  if (D==2){
    sigma = diag(D)
    sigma[sigma==0] = 0.5
  }else{
    sigma = randcorr(D)
  }
  
  if (itemnum == 4 & testlet == 4){
    sigma = matrix(0,4,4)
    sigma[1,2] = 0.3
    sigma[1,3] = 0.1
    sigma[1,4] = 0.3
    sigma[2,3] = 0.3
    sigma[2,4] = 0.1
    sigma[3,4] = 0.3
    sigma = sigma + t(sigma)
    diag(sigma) <- 1
  }
  
  
  set.seed(rand.seed)
  
  theta <- rmvnorm(n=people, mean=rep(0,D), sigma=sigma)
  
  
  

  
  # simulate responses
  MU = rep(0,itemnum)
  
  # item_index = matrix(NA,testlet,2)
  # item_index[1,1] = 1
  # item_index[1,2] = itemnum
  # for (k in 2:testlet){
  #   item_index[k,1] = item_index[k-1,1] + itemnum
  #   item_index[k,2] = item_index[k-1,2] + itemnum
  # }
  
  item_index = t(matrix(1:total_items,D,testlet))
  
  
  if (simulation == TRUE & is.null(response)){
    if ( (use_real_data_pars == 0 | model == 0 | model == 4) & use_real_data_pars <= 0){
      
      res = vector("list",testlet)
      for (k in 1:testlet){
        
        
        # # model identification (I constrain variance=1 of latent traits)
        # if (k == 1) alpha[[k]][] = 1 
        
        # no correlations for errors
        # error[[k]] <- rmvnorm(n=people, mean=rep(0,D), sigma=diag(D))
        
        # simulate response
        res[[k]] = matrix(NA,people,itemnum)
        for(n in 1:people){
          for (i in 1:itemnum){
            res[[k]][n,i] = delta[[k]][i] + alpha[[k]][i,]%*%theta[n,] + error[[k]][n,i]
          }
        }

      }
      y = do.call(cbind,res)
      
      if (model == 0){
        # inverse logit function: Transforming (-Inf,Inf) to original scale (0,1)
        response = plogis(y) # formula is 1/(1 + exp(-u))
      }else if (model == 4){
        response = y
      }
      
    }else if (use_real_data_pars == 1 | use_real_data_pars == 2 | use_real_data_pars == 3 
               | model == 1 | model == 2 | model == 3){
      
      response <- matrix(NA,people,itemnum*testlet) # testlet * itemnum
      for(n in 1:people){
        cnt = 1:itemnum
        for (k in 1:testlet){
          
          M = N = rep(NA,itemnum)
          for(i in 1:itemnum){
            M[i] <- exp(((alpha[[k]][i,]%*%theta[n,] - delta[[k]][i])+tau[[k]][i])/2)
            N[i] <- exp((-(alpha[[k]][i,]%*%theta[n,] - delta[[k]][i])+tau[[k]][i])/2)
          }
          
          z = rmvnorm(n=1, mean=MU, sigma=sigma_item_cor[[k]])
          response[n,cnt] = qbeta(pnorm(z, 0, sqrt(diag(sigma_item_cor[[k]]))), M, N)
          
          while (any(response[n,cnt]>=1)){
            z = rmvnorm(n=1, mean=MU, sigma=sigma_item_cor[[k]])
            response[n,cnt] = qbeta(pnorm(z, 0, sqrt(diag(sigma_item_cor[[k]]))), M, N)
            print("response >= 1.")
          }
          cnt = cnt + itemnum
        }
      }
      
    }

  }

  
  
  
  
  # saveRDS(response, file = "y01.rds")
  # response = readRDS(file = "y01.rds")
  
  
  
if (return_sim==TRUE){
  if (model == 0 | model == 4) return(list(response=response
                              ,alpha=alpha
                              ,delta=delta
                              ,sigma_item_cor=sigma_item_cor
                              ,error=error
                              ,theta=theta
                              ,sigma=sigma
                              ))
  if (model == 1 | model == 2 | model == 3 | model == 4){
    return(list(response=response
                ,alpha=alpha
                ,delta=delta
                ,tau=tau
                ,theta=theta
                ,sigma=sigma
                ,sigma_item_cor=sigma_item_cor
                ))
  }
}
  
  
  
  
  
  
  
  # nimble code
  require("nimble")
  require("evaluate")
  

  
  
  # for posterior predictive sampling
  ppSamplerNF <- nimbleFunction(
    setup = function(model, mcmc) {
      dataNodes <- model$getNodeNames(dataOnly = TRUE)
      parentNodes <- model$getParents(dataNodes, stochOnly = TRUE)
      # cat("Stochastic parents of data are:", paste(parentNodes, collapse = ','), ".\n")
      simNodes <- model$getDependencies(parentNodes, self = FALSE)
      vars <- mcmc$mvSamples$getVarNames()  # need ordering of variables in mvSamples / samples matrix
      # cat("Using posterior samples of:", paste(vars, collapse = ','), ".\n")
      n <- length(model$expandNodeNames(dataNodes, returnScalarComponents = TRUE))
    },
    run = function(samples = double(2)) {
      nSamp <- dim(samples)[1]
      ppSamples <- matrix(nrow = nSamp, ncol = n)
      for(i in 1:nSamp) {
        values(model, vars) <<- samples[i, ]
        model$simulate(simNodes, includeData = TRUE)
        ppSamples[i, ] <- values(model, dataNodes)
      }
      returnType(double(2))
      return(ppSamples)
    })
  
  
  
  
  if (model == 0){ # logit-normal model
    
    text = paste('

code <- nimbleCode({

  L[1:D,1:D] ~ dlkj_corr_cholesky(eta=eta, p=D)
  corr[1:D,1:D] <- t(L[1:D,1:D])%*%L[1:D,1:D] # correlation matrix
  
  for (k in 1:testlet){
    Ustar[1:D,1:D,k] ~ dlkj_corr_cholesky(eta=eta, p=D)
    U[1:D,1:D,k] <- uppertri_mult_diag(Ustar[1:D,1:D,k], sds[item_index[k,1]:item_index[k,D]])
    sigma_item_cor[1:D,1:D,k] <- t(U[1:D,1:D,k])%*%U[1:D,1:D,k] # covariance matrix
  }
  

  for(n in 1:people){
  
    for (k in 1:total_items){
      MU[n,k] <- delta[k] + inprod(alpha[k,1:D],theta[n,1:D])
    }
    for (k in 1:testlet){
      response[n,item_index[k,1]:item_index[k,itemnum]] ~ dlogitnormal(MU[n,item_index[k,1]:item_index[k,itemnum]], cholesky = U[1:D,1:D,k], prec_param = 0)
    }
    theta[n, 1:D] ~ dmnorm(mu[1:D], cholesky = L[1:D, 1:D], prec_param = 0)
    
  }

  for (k in 1:total_items){
    delta[k] ~  dnorm(0,0.0001)
    E2[k] ~ dnorm(0,1)
    D2[k] ~ dnorm(0,0.0016)
    sds[k] <- pow(E2[k]/D2[k],2) # Half-Cauchy prior
  }
  
  for (k in 1:index_alpha_length){
    E1[k] ~ dnorm(0,1)
    D1[k] ~ dnorm(0,0.0016)
    alpha[index_alpha[k,1],index_alpha[k,2]] <- pow(E1[k]/D1[k],2) # Half-Cauchy prior
  }
  for (k in 1:index_zero_length){
    alpha[index_zero[k,1],index_zero[k,2]] <- 0
  }

})

')
    
  }else if (model == 1){  # multidimensional beta model
    
    text = paste('

code <- nimbleCode({

  L[1:D,1:D] ~ dlkj_corr_cholesky(eta=eta, p=D)

  for(n in 1:people){
    for (k in 1:total_items){
      M[n,k] <- exp(((inprod(alpha[k,1:D],theta[n,1:D])-delta[k])+tau[k])/2)
      N[n,k] <- exp((-(inprod(alpha[k,1:D],theta[n,1:D])-delta[k])+tau[k])/2)
      response[n,k] ~ dbeta(M[n,k],N[n,k])
    }
    theta[n, 1:D] ~ dmnorm(mu[1:D], cholesky = L[1:D, 1:D], prec_param = 0)
  }
  
  corr[1:D,1:D] <- t(L[1:D,1:D])%*%L[1:D,1:D] # correlation matrix

  for (k in 1:total_items){
    delta[k] ~ dnorm(0,0.0001)
    tau[k] ~  dnorm(0,0.0001)
  }
  for (k in 1:index_alpha_length){
    E1[k] ~ dnorm(0,1)
    D1[k] ~ dnorm(0,0.0016)
    alpha[index_alpha[k,1],index_alpha[k,2]] <- pow(E1[k]/D1[k],2) # Half-Cauchy prior
  }
  for (k in 1:index_zero_length){
    alpha[index_zero[k,1],index_zero[k,2]] <- 0
  }
})

')
    
    
    
    
  }else if (model == 2){ # multidimensional beta model with item correlations
    
    text = paste('

code <- nimbleCode({

  for(n in 1:people){
    for (k in 1:total_items){
      M[n,k] <- exp(((inprod(alpha[k,1:D],theta[n,1:D])-delta[k])+tau[k])/2)
      N[n,k] <- exp((-(inprod(alpha[k,1:D],theta[n,1:D])-delta[k])+tau[k])/2)
    }
    for (k in 1:testlet){
      response[n,item_index[k,1]:item_index[k,itemnum]] ~ dmbeta(a = M[n,item_index[k,1]:item_index[k,itemnum]], b = N[n,item_index[k,1]:item_index[k,itemnum]], cholesky = L22[1:D,1:D,k], prec_param = 0)
    }
    theta[n, 1:D] ~ dmnorm(mu[1:D], cholesky = L[1:D, 1:D], prec_param = 0)
  }
  
  L[1:D,1:D] ~ dlkj_corr_cholesky(eta=eta, p=D)
  corr[1:D,1:D] <- t(L[1:D,1:D])%*%L[1:D,1:D] # correlation matrix
  
  for (k in 1:testlet){
    L22[1:D,1:D,k] ~ dlkj_corr_cholesky(eta=eta, p=D)
    sigma_item_cor[1:D,1:D,k] <- t(L22[1:D,1:D,k])%*%L22[1:D,1:D,k] # correlation matrix
  }

  for (k in 1:total_items){
    delta[k] ~ dnorm(0,0.0001)
    tau[k] ~  dnorm(0,0.0001)
  }
  for (k in 1:index_alpha_length){
    E1[k] ~ dnorm(0,1)
    D1[k] ~ dnorm(0,0.0016)
    alpha[index_alpha[k,1],index_alpha[k,2]] <- pow(E1[k]/D1[k],2) # Half-Cauchy prior
  }
  for (k in 1:index_zero_length){
    alpha[index_zero[k,1],index_zero[k,2]] <- 0
  }
})

')
    
    
    
    
  }else if (model == 3){ # multidimensional beta model without item correlations
    
    
    
    text = paste('

code <- nimbleCode({

  for(n in 1:people){
    for (k in 1:total_items){
      M[n,k] <- exp(((inprod(alpha[k,1:D],theta[n,1:D])-delta[k])+tau[k])/2)
      N[n,k] <- exp((-(inprod(alpha[k,1:D],theta[n,1:D])-delta[k])+tau[k])/2)
    }
    for (k in 1:testlet){
      response[n,item_index[k,1]:item_index[k,itemnum]] ~ dmbeta(a = M[n,item_index[k,1]:item_index[k,itemnum]], b = N[n,item_index[k,1]:item_index[k,itemnum]], cholesky = L22[1:D,1:D,k], prec_param = 0)
    }
    theta[n, 1:D] ~ dmnorm(mu[1:D], cholesky = L[1:D, 1:D], prec_param = 0)
  }
  
  L[1:D,1:D] ~ dlkj_corr_cholesky(eta=eta, p=D)
  corr[1:D,1:D] <- t(L[1:D,1:D])%*%L[1:D,1:D] # correlation matrix
  
  # for (k in 1:testlet){
  #   # L22[1:D,1:D,k] ~ dlkj_corr_cholesky(eta=eta, p=D)
  #   sigma_item_cor[1:D,1:D,k] <- t(L22[1:D,1:D,k])%*%L22[1:D,1:D,k] # correlation matrix
  # }

  for (k in 1:total_items){
    delta[k] ~ dnorm(0,0.0001)
    tau[k] ~  dnorm(0,0.0001)
  }
  for (k in 1:index_alpha_length){
    E1[k] ~ dnorm(0,1)
    D1[k] ~ dnorm(0,0.0016)
    alpha[index_alpha[k,1],index_alpha[k,2]] <- pow(E1[k]/D1[k],2) # Half-Cauchy prior
  }
  for (k in 1:index_zero_length){
    alpha[index_zero[k,1],index_zero[k,2]] <- 0
  }
})

')
    
    
  }else if (model == 4){  # correlated traits–correlated uniqueness model (CTCUM; Marsh, 1989)
    
    text = paste('

code <- nimbleCode({

  L[1:D,1:D] ~ dlkj_corr_cholesky(eta=eta, p=D)
  
  for (k in 1:testlet){
    Ustar[1:D,1:D,k] ~ dlkj_corr_cholesky(eta=eta, p=D)
    U[1:D,1:D,k] <- uppertri_mult_diag(Ustar[1:D,1:D,k], sds[item_index[k,1]:item_index[k,D]])
    sigma_item_cor[1:D,1:D,k] <- t(U[1:D,1:D,k])%*%U[1:D,1:D,k] # covariance matrix
  }
  

  for(n in 1:people){
    for (k in 1:testlet){
      for (i in 1:itemnum){
        MU[n,i,k] <- delta[item_index[k,i]] + alpha[item_index[k,i]]*theta[n,i]
      }
      response[n,item_index[k,1]:item_index[k,D]] ~ dmnorm(MU[n,1:D,k], cholesky = U[1:D,1:D,k], prec_param = 0)
    }
    theta[n, 1:D] ~ dmnorm(mu[1:D], cholesky = L[1:D, 1:D], prec_param = 0)
  }
  
  corr[1:D,1:D] <- t(L[1:D,1:D])%*%L[1:D,1:D] # correlation matrix

  for (k in 1:total_items){
    E1[k] ~ dnorm(0,1)
    D1[k] ~ dnorm(0,0.0016)
    alpha[k] <- pow(E1[k]/D1[k],2) # Half-Cauchy prior
    delta[k] ~  dnorm(0,0.0001)
    E2[k] ~ dnorm(0,1)
    D2[k] ~ dnorm(0,0.0016)
    sds[k] <- pow(E2[k]/D2[k],2) # Half-Cauchy prior
  }

})

')
    
  }
  

  
  
  # nimble settings
  ex = evaluate(text)
  
  if (model == 0 | model == 4){ # logit-normal model
    
    constants <- list(people=people,itemnum=itemnum
                      ,testlet=testlet,D=D,mu=mu
                      ,eta=1
                      ,total_items=total_items
                      ,item_index=item_index
                      ,index_alpha=index_alpha
                      ,index_zero=index_zero
                      ,index_alpha_length=index_alpha_length
                      ,index_zero_length=index_zero_length)
    
    
    
    data <- list(response=response)
    
    
    # nimble is sensitive to initial values
    inits <- list(E1=rep(0.5,total_items),D1=rep(0.5,total_items)
                  ,E2=rep(0.5,total_items),D2=rep(0.5,total_items)
                  ,delta=rep(0,total_items)
                  ,L=diag(D)
                  ,Ustar=Ustar
    )
    
    monitors = c("alpha","delta","corr","sigma_item_cor")
    
  }else{
    
    constants <- list(people=people,itemnum=itemnum
                      ,testlet=testlet,D=D,mu=mu
                      ,eta=eta
                      ,total_items=total_items
                      ,index_alpha=index_alpha
                      ,index_zero=index_zero
                      ,index_alpha_length=index_alpha_length
                      ,index_zero_length=index_zero_length
                      ,item_index=item_index)
    
    
    
    data <- list(response=response)
    
    
    # nimble is sensitive to initial values
    inits <- list(E1=rep(0.5,total_items),D1=rep(0.5,total_items) #alpha=alpha_ini,
                  ,E2=rep(0.5,total_items),D2=rep(0.5,total_items)
                  ,delta=rep(0,total_items)
                  ,tau=rep(0,total_items)
                  ,L=diag(D)
                  ,L22=L22
    ) # ,theta=matrix(0,people,D)
    
    
    monitors = c("alpha","delta","tau","corr")
    
    if (model == 2 | model == 4){
      monitors = c(monitors, "sigma_item_cor" )
    }
  }
  



  
  if (simulation == FALSE | est_theta == TRUE | sim_data == TRUE){
    monitors_add = c()
    monitors_add = c(monitors_add, "theta" )
    
    if ( (model == 0) & simulation == FALSE){
      monitors_add = c(monitors_add, "E1", "D1", "E2", "D2", "Ustar")
    }
    if ( (model == 1 | model == 2 | model == 3) & simulation == FALSE){
      monitors_add = c(monitors_add, "E1", "D1")
      # monitors_add = c(monitors_add, "M", "N") # for residual calculation
    }
    if (model == 2 & simulation == FALSE){
      monitors_add = c(monitors_add, "L22")
    }
    if ( (model == 4) & simulation == FALSE){
      monitors_add = c(monitors_add, "E1", "D1", "E2", "D2")
    }
    
    if (simulation == FALSE){
      monitors_add = c(monitors_add, "logProb_response")
    }
    
    monitors_vital = monitors
    monitors = c(monitors, monitors_add)
  }
  
  enableWAIC = WAIC = TRUE
  if (simulation == FALSE){
    enableWAIC = WAIC =FALSE
  }
  if (waic_nimble==TRUE){
    enableWAIC = WAIC = TRUE
  }
  
  
  MODEL <- nimbleModel(code = code, name = "MODEL", constants = constants,data = data, inits = inits)

  # generating C++ code
  cmodel <- compileNimble(MODEL)

  conf <- configureMCMC(MODEL, monitors = monitors
                        ,useConjugacy=TRUE ,print = FALSE
                        ,enableWAIC=enableWAIC)

  # conf$printSamplers("delta[1]") # what samplers used

  modelMCMC <- buildMCMC(conf)
  cmodelMCMC <- compileNimble(modelMCMC, project = MODEL ,resetFunctions = TRUE)

  niter = burnin + keep

  

  # Check if Geweke’s convergence diagnostic is satisfied. If not, continue the MCMC chains. 
  library("coda")
  sig1_items = sig2_items = TRUE
  if (nchains==1) sig2_items = FALSE
  if (nchains>1) sig1_items = FALSE
  cnt = 1
  
  if (nchains>1 & simulation==TRUE & mcmc_length_add != 0) stop("How to use cmodelMCMC$run for multiple chains!")
  
  while (any(sig1_items) | any(sig2_items)){
    
    if (cnt >=2) print(paste("Continue MCMC chains it left out. ",cnt," times.",sep=""))
    
    if (cnt == 1){
      times <- system.time(samples <- runMCMC(cmodelMCMC, niter = niter, nburnin = burnin, nchains=nchains
                                              ,samplesAsCodaMCMC=TRUE
                                              ,samples=TRUE
                                              ,summary=TRUE
                                              ,WAIC = WAIC
                                              ,setSeed = rand.seed))
      if ( mcmc_length_add != 0) ma <- as.matrix(cmodelMCMC$mvSamples)
    }else{
      cmodelMCMC$run(mcmc_length_add, reset = FALSE, resetMV = TRUE, time = FALSE, progressBar = FALSE, nburnin = 0, thin = -1, thin2 = -1, resetWAIC = FALSE) # run for more samples
      ma <- rbind(ma, as.matrix(cmodelMCMC$mvSamples))
      end = keep + mcmc_length_add
      start = (end - keep) + 1
      ma = ma[start:end, ]
      samples$samples = ma
      samples$summary <- samplesSummary(samples$samples) # update summary
      samples$WAIC = cmodelMCMC$getWAIC() # update WAIC (NEED ATTENTION)
    }

    
    sam = samples$samples
    if (nchains == 1) sam = list(sam)
    samples_items <- samples_theta <- vector("list", length(sam))
    for (k in 1:length(sam)){
      for (j in 1:length(monitors_vital)){
        ind <- grepl(monitors_vital[j], colnames(sam[[k]]))
        samples_items[[k]] <- cbind(samples_items[[k]], sam[[k]][,ind])
      }
      samples_items[[k]] <- mcmc(samples_items[[k]])
      
      ind <- grepl("theta", colnames(sam[[k]]))
      samples_theta[[k]] <- mcmc(sam[[k]][,ind])
    }
    
    # obtain gelman.diag
    if (nchains>=2){
      gelman_item <- gelman.diag(samples_items, confidence = 0.95, transform=FALSE, autoburnin=FALSE, multivariate=FALSE)
      # gelman_theta <- gelman.diag(samples_theta, confidence = 0.95, transform=FALSE, autoburnin=FALSE, multivariate=FALSE)
      
      value_items = gelman_item$psrf[,2] # PSRF Upper C.I.
      sig2_items = value_items[!is.nan(value_items)] > 1.1
      nonconverged_items = value_items > 1.1
      
      value_theta = gelman_theta$psrf[,2] # PSRF Upper C.I.
      sig2_theta = value_theta[!is.nan(value_theta)] > 1.1
      nonconverged_theta = value_theta > 1.1
    }
    
    # Geweke’s convergence diagnostic
    if (nchains == 1){
      geweke_items = geweke.diag(samples_items, frac1=0.1, frac2=0.5)
      pvalue_items = pnorm(abs(geweke_items$z),lower.tail=FALSE)*2 # If pvalue < .05, it does not converge yet!
      sig1_items = pvalue_items[!is.nan(pvalue_items)] < .05
      nonconverged_items = pvalue_items < .05
      
      geweke_theta = geweke.diag(samples_theta, frac1=0.1, frac2=0.5)
      pvalue_theta = pnorm(abs(geweke_theta$z),lower.tail=FALSE)*2 # If pvalue < .05, it does not converge yet!
      sig1_theta = pvalue_theta[!is.nan(pvalue_theta)] < .05
      nonconverged_theta = pvalue_theta < .05
    }
    
    # Heidelberger and Welch’s convergence diagnostic
    if (nchains == 1){
      HW_items = heidel.diag(samples_items[[1]], eps=0.1, pvalue=0.05)
      HW_theta = heidel.diag(samples_theta[[1]], eps=0.1, pvalue=0.05)
    }

    # If analyzing real data and mcmc_length_add==0, we run MCMC once. 
    if (simulation == FALSE) break
    if (mcmc_length_add==0) break
    
    cnt = cnt + 1
  }
  
  samples$mcmc_length_add_times = cnt
  samples$nonconverged_items = nonconverged_items
  samples$nonconverged_theta = nonconverged_theta
  
  samples$HW_items = HW_items
  samples$HW_theta = HW_theta
  
  # combine chains for later use
  if (nchains>=2){
    samples$summary = samples$summary$all.chains
    samples$samples = do.call(rbind,samples$samples)
  }
  
  # times <- system.time(samples <- nimbleMCMC(code = code, #model = MODEL,
  #                                   constants = constants,
  #                                   data = data, inits = inits,
  #                                   nchains = 3,
  #                                   nburnin = 10000, niter = 40000,
  #                                   WAIC = TRUE,
  #                                   monitors = monitors,
  #                                   summary = TRUE,
  #                                   setSeed = TRUE,
  #                                   samplesAsCodaMCMC = TRUE))
  
  
  # cmodelMCMC$getWAIC()
  # cmodelMCMC$getWAICdetails()
  
  if (get_item_fit == TRUE & model==0){
    print("Item fit for logit-normal model is not built because calculating mean adn variance are difficult.")
  }
  
  if (get_item_fit == TRUE & (model==0 | model == 1 | model == 2)){
    
    # times_mn <- system.time()
    
    mn <- beta_M_N(response, model, item_index
                   , samples$samples
                   , keep, people
                   , itemnum, testlet, total_items)
    
    # Standardized generalized dimensionality discrepancy measure (SGDDM)
    samples$cov_obs = mn$cov_obs
    samples$cov_rep = mn$cov_rep
    if (model == 1 | model == 2){
      samples$cov_ppp <- apply(mn$cov_rep > mn$cov_obs, c(2,3), mean)
      
      samples$sgddm_obs = mn$sgddm_obs
      samples$sgddm_rep = mn$sgddm_rep
    }
    
    # # skewness and ...
    # samples$st = mn$st
    
    
    samples$prob = mn$prob
    samples$r_Cox_Snell = mn$r_Cox_Snell
    samples$eap_r_Cox_Snell = mn$eap_r_Cox_Snell
    samples$pred_lower_yi = mn$pred_lower_yi
    
    samples$eap_sigma_item_cor_4D = mn$eap_sigma_item_cor_4D
    samples$eap_alpha = mn$eap_alpha
    samples$eap_theta = mn$eap_theta
    samples$eap_delta = mn$eap_delta
    samples$eap_tau = mn$eap_tau
    
    if (model == 1 | model == 2){
      z2_obs = mn$nu_obs2/mn$var
      z2_rep = mn$nu_rep2/mn$var
      
      # outfit: item fit
      r_obs_outfit_item <- apply(z2_obs,c(1,3),mean,na.rm=TRUE)
      r_rep_outfit_item <- apply(z2_rep,c(1,3),mean,na.rm=TRUE)
      ppp_outfit_item = apply(r_rep_outfit_item > r_obs_outfit_item, 2, mean)
      
      # outfit: person fit
      r_obs_outfit_person <- apply(z2_obs,c(1,2),mean,na.rm=TRUE)
      r_rep_outfit_person <- apply(z2_rep,c(1,2),mean,na.rm=TRUE)
      ppp_outfit_person = apply(r_rep_outfit_person > r_obs_outfit_person, 2, mean)
      
      
      # infit: item fit
      r_obs_infit_item <- apply(mn$nu_obs2,c(1,3),mean,na.rm=TRUE) / apply(mn$var,c(1,3),mean,na.rm=TRUE)
      r_rep_infit_item <- apply(mn$nu_rep2,c(1,3),mean,na.rm=TRUE) / apply(mn$var,c(1,3),mean,na.rm=TRUE)
      ppp_infit_item = apply(r_rep_infit_item > r_obs_infit_item, 2, mean)
      
      # infit: person fit
      r_obs_infit_person <- apply(mn$nu_obs2,c(1,2),mean,na.rm=TRUE) / apply(mn$var,c(1,2),mean,na.rm=TRUE)
      r_rep_infit_person <- apply(mn$nu_rep2,c(1,2),mean,na.rm=TRUE) / apply(mn$var,c(1,2),mean,na.rm=TRUE)
      ppp_infit_person = apply(r_rep_infit_person > r_obs_infit_person, 2, mean)
    }
    
    rm(mn)
    gc(reset=TRUE)
  }
  
  # obtain reliability of latent traits:
  reliability = NULL
  if (simulation == FALSE | est_theta == TRUE){
    ind_theta <- grepl('theta', rownames(samples$summary))
    ind_logProb_theta <- grepl('logProb_theta', rownames(samples$summary))
    ind_theta[ind_logProb_theta] = FALSE # ignore string "logProb_theta"
    theta_est = samples$summary[ind_theta,"Mean"]
    theta_sd = samples$summary[ind_theta,"St.Dev."]
    dim(theta_est) <- c(people, D)
    dim(theta_sd) <- c(people, D)
    reliability = apply(theta_est,2,var) / (apply(theta_est,2,var) + apply(theta_sd^2,2,mean))
  }
  
  
  if (simulation == FALSE){
    # keep33 = keep
    # logProb = as.matrix(cmodelMCMC$mvSamples) # for only one chain
    
    keep33 = keep * nchains
    logProb = samples$samples # for multiple chains
    ind <- grepl('logProb_response', colnames(logProb))
    logProb_response <- logProb[,ind]
    if (model == 0 | model == 2 | model == 3 | model == 4){
      # because of mvn, the logProb will have only one value for a vector of responses
      dim(logProb_response) <- c(keep33,people,item_index[testlet,1])
      logProb_response = logProb_response[,,item_index[,1]]
    }else if (model == 1){
      dim(logProb_response) <- c(keep33,people,total_items)
    }
    
    # WAIC : see https://stats.stackexchange.com/questions/331928/dic-waic-in-jags
    Prob_response = exp(logProb_response)
    Prob_response_mean = apply(Prob_response,c(2,3),mean)
    lppd = sum(log(Prob_response_mean)) # log pointwise predictive density
    # pwaic2
    pwaic2 = sum(apply(logProb_response,c(2,3),var))
    # WAIC 2
    waic2 = -2*(lppd-pwaic2)

    samples$lppd = lppd
    samples$pwaic2 = pwaic2
    samples$waic2 = waic2
    

    
  }
  
  

  
  

  options(max.print=1000000)

  
  if (simulation == TRUE){
    samples$samples = NULL
    if (get_item_fit==TRUE) samples$response = response
  }
  samples$model_is = model_is
  samples$times = times
  samples$itemnum = itemnum
  samples$testlet = testlet
  samples$total_items = total_items
  samples$people = people
  samples$burnin = burnin
  samples$keep = keep
  samples$item_index = item_index
  samples$nchains = nchains
  
  if (simulation == TRUE){
    samples$alpha = alpha
    samples$delta = delta
    if (model == 1 | model == 2 | model == 3) samples$tau = tau
    samples$corr = sigma
    samples$sigma_item_cor = sigma_item_cor
    samples$theta = theta
  }
  
  if ( (simulation == FALSE) | est_theta == TRUE | sim_data == TRUE){
    samples$monitors_vital = monitors_vital
    
    monitors = monitors[ -which(monitors == "logProb_response")]
    samples$monitors_add = monitors_add[ -which(monitors == "logProb_response")]
    if (simulation == FALSE) {
      samples$logProb_response = logProb_response
    }
  }
  

  # extract all summary of parameters of interest
  samples$est = list()
  for (k in 1:length(samples$monitors_vital)){
    name = samples$monitors_vital[k]
    name_all = rownames(samples$summary)
    index <- grepl(name, name_all)
    samples$est[[name]] = samples$summary[index,]
  }
  if (est_theta==TRUE){
    samples$est$theta_est = theta_est
    samples$est$theta_se = theta_sd
  }


  # test reliability
  samples$reliability = reliability
  # samples$geweke = geweke
  
  
  samples$code = code
  samples$constants = constants
  if (get_item_fit==TRUE) samples$data = data
  samples$inits = inits
  samples$monitors = monitors
  
  if (get_item_fit == TRUE & (model == 1 | model == 2)){
    # samples$r_obs_outfit_item = r_obs_outfit_item
    # samples$r_rep_outfit_item = r_rep_outfit_item
    samples$ppp_outfit_item = ppp_outfit_item
    
    # samples$r_obs_outfit_person = r_obs_outfit_person
    # samples$r_rep_outfit_person = r_rep_outfit_person
    samples$ppp_outfit_person = ppp_outfit_person
    
    # samples$r_obs_infit_item = r_obs_infit_item
    # samples$r_rep_infit_item = r_rep_infit_item
    samples$ppp_infit_item = ppp_infit_item
    
    # samples$r_obs_infit_person = r_obs_infit_person
    # samples$r_rep_infit_person = r_rep_infit_person
    samples$ppp_infit_person = ppp_infit_person
  }
  
  if (nchains>1){
    samples$samples_items = samples_items
    samples$samples_theta = samples_theta
  }
  
  return(samples)
  
}


