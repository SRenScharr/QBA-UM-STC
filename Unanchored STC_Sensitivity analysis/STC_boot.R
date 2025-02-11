#A function for unanchored STC analysis for a binary outcome
#Note that the default settings have 0 continuous variable for both measured and unmeasured covariates
#Note this script need to be named as “STC_boot.R”

STC_boot <- function(IPD_data,
                     AgD_data,
                     measured_cov, #the measured covariates
                     unmeasured_cov=NULL,
                     con_ind=0, #a number to indicate the location of measured continuous covariate, take value 0 if null
                     un_con_ind=0, #a number to indicate the location of unmeasured continuous covariate, take value 0 if null
                     un_con_mean=0,
                     un_con_sd=0,
                     un_bin_mean=0,
                     n_b=10^4, #number of bootstrap
                     n_me=10^4){ #number of covariates to sample for measured covariates
  
  #Assuming the AgD study had treatment A
  #Summary statistics for arm A
  r_A <- AgD_data$number.of.events[1]
  n_A <- AgD_data$number.of.patients[1]
  d_A <- log(r_A/n_A/(1-r_A/n_A))
  var_A <- n_A/(r_A*(n_A-r_A))
  
  #Construct the copula structure for measured_cov for simulating covariates 
  if(con_ind[1]==0) { #if there is no continuous covariates, all variables are binary
    cov_con <- NULL
    cov_bin <- measured_cov
  }else{ #if there is at least one continuous covariates
    cov_con <- measured_cov[con_ind]
    cov_bin <- measured_cov[-con_ind] 
  }
  cov_name_int <- c(cov_con,cov_bin) #reorder covariates with continuous first
  n_con <- length(cov_con)
  n_bin <- length(cov_bin)
  margin_type <- c(rep("norm", n_con),rep("binom", n_bin))
  con_list <- vector("list", length =  n_con) #define an empty list
  bin_list <- vector("list", length =  n_bin) #define an empty list
  
  
  #Construct the copula structure for unmeasured_cov for simulating covariates 
  if(un_con_ind[1]==0) {
    un_cov_con <- NULL
    un_cov_bin <- unmeasured_cov
  }else{
    un_cov_con <- unmeasured_cov[un_con_ind]
    un_cov_bin <- unmeasured_cov[-un_con_ind] 
  }
  un_cov_name_int <- c(un_cov_con,un_cov_bin) #reorder covariates with continuous first
  un_n_con <- length(un_cov_con)
  un_n_bin <- length(un_cov_bin)
  un_margin_type <- c(rep("norm",un_n_con),rep("binom", un_n_bin))
  un_con_list <- vector("list", length=un_n_con) #define an empty list
  un_bin_list <- vector("list", length=un_n_bin) #define an empty list
  
  #The name of the full covariate list (including both meansured and unmeasured)
  full_cov <- c(cov_name_int,un_cov_name_int)
  
  #Assuming the IPD study had treatment B
  #Bootstrap data for arm B using the IPD data
  boot_fun <- function(i) {
    index <- sample(nrow(IPD_data), nrow(IPD_data), replace = TRUE)
    boot_IPD <- IPD_data[index, ]
    ipd_X<- boot_IPD[,full_cov]
    ipd_cor <- cor(ipd_X)
    ipd_copula <-
      copula::normalCopula(copula::P2p(ipd_cor),
                           dim = ncol(ipd_cor),
                           dispstr = "un") #construct a normal copula for n_X random variables with correlation structure as ipd_cor
    
    if (n_con>0){
      for (ic in 1:n_con){
        con_list[[ic]] <- list(mean = AgD_data[cov_con[ic]][1,], sd = AgD_data[cov_con[ic]][2,])
      } }
    
    if (n_bin>0){
      for (ib in 1:n_bin){
        bin_list[[ib]] <- list(size=1, prob = AgD_data[cov_bin[ib]][1,])
      } }
    if (un_n_con>0){
      for (ic in 1:un_n_con){
        un_con_list[[ic]] <- list(mean = un_con_mean[ic], sd = un_con_sd[ic])
      } }
    
    if (un_n_bin>0){
      for (ib in 1:un_n_bin){
        un_bin_list[[ib]] <- list(size=1, prob = un_bin_mean[ib])
      } }
    
    Mvd <-
      copula::mvdc(
        copula = ipd_copula,
        margins = c(margin_type,un_margin_type),
        paramMargins = c(con_list,bin_list,un_con_list,un_bin_list)
      )
    
    AgD_cov <- copula::rMvdc(n_me, Mvd) # simulated covariates for the AgD population
    AgD_cov <- data.frame(AgD_cov)
    colnames(AgD_cov) <- full_cov
    
    #Obtain a regression model using the IPD study
    m_STC1 <-
      glm(as.formula(paste("y ~ ",
                           paste(
                             full_cov, collapse = "+"
                           ))), data = boot_IPD, family = binomial)
    
    #Predict outcome for the AgD population
    p_STC1 <- predict(m_STC1,
                      newdata = AgD_cov,
                      type = "response")
    p_B <- mean(p_STC1)
    
    log(p_B/ (1 - p_B))
  }
  
  boot.results <- sapply(1:n_b, boot_fun)
  
  log_OR <- mean(boot.results)-d_A
  log_OR_l <- log_OR-qnorm(0.975)*sqrt(var(boot.results)+var_A)
  log_OR_u <- log_OR+qnorm(0.975)*sqrt(var(boot.results)+var_A)
  
  #return point estimate and 95%CI
  c(exp(log_OR),exp(log_OR_l),exp(log_OR_u))
}
