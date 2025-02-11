#Simulate the IPD and AgD study to be used as an example
#Step1: Inputs for the IPD study
#Step2: Inputs for the AgD study
#Step3: Simulate the IPD study
#Step4: Simulate the AgD study
#Step5: Data used in the STC

#Load library
library(dplyr)
library(copula)

#Set working directory
#path <- ""
#setwd(path)

#Set seeds
set.seed(12345)

################################################################################################
#STEP1: Specify inputs for simulating the IPD study
################################################################################################
N_IPD <- 300 #sample size 
covariates_IPD <- c("X1", "X2","X3","X4") #name of the covariates: two continuous and two binary 
p_IPD<- length(covariates_IPD)
cor_input_IPD<- matrix(0,nrow=p_IPD,ncol=p_IPD) #correlation structure
diag(cor_input_IPD) <- 1
m_IPD <- matrix(c(23,1, 75,2, 0.3,0, 0.75,0),nrow=2) #summary statistics for the covariates
rownames(m_IPD) <- c("mean","sd")
b <- c(0.05,0.008,-0.002,-0.05,1) #coefficients (intercept, X1, X2, X3, X4)
b_trt <- -0.5

################################################################################################
#STEP2: Specify inputs for simulating the AgD study 
################################################################################################
N_AgD <- 300 #sample size
covariates_AgD <- c("X1", "X2","X3","X4") #name of the covariates: two continuous and two binary 
p_AgD<- length(covariates_AgD)
cor_input_AgD<- matrix(0,nrow=p_AgD,ncol=p_AgD) #correlation structure
diag(cor_input_AgD) <- 1
m_AgD <- matrix(c(24,1, 76,2, 0.2,0, 0.75,0),nrow=2) #summary statistics for the covariates
rownames(m_AgD) <- c("mean","sd")


################################################################################################
#STEP3: Use NORTA/Gaussian copula to simulate covariates for IPD study
################################################################################################
myMvd_IPD <-
  mvdc(
    copula = normalCopula(copula::P2p(cor_input_IPD), dim = p_IPD, dispstr = "un"),
    margins = c("norm","norm","binom","binom"),
    paramMargins = list(list(mean = m_IPD[1,1], sd = m_IPD[2,1]),
                        list(mean = m_IPD[1,2], sd = m_IPD[2,2]),
                        list(size = 1, prob = m_IPD[1,3]),
                        list(size = 1, prob = m_IPD[1,4]))
  )
mycov_IPD <-rMvdc(N_IPD, myMvd_IPD)

#Simulate binary outcome 
bX_IPD <- b[1]+b[-1]%*%t(mycov_IPD)+b_trt
y_IPD <- rbinom(N_IPD,1,1/(1+exp(-bX_IPD)))

#Simulated IPD for the IPD study
IPD <- data.frame(cbind(y_IPD,mycov_IPD))
colnames(IPD) <- c("y",covariates_IPD)

################################################################################################
#STEP4: Use NORTA/Gaussian copula to simulate covariates for AgD study
################################################################################################
myMvd_AgD <-
  mvdc(
    copula = normalCopula(copula::P2p(cor_input_AgD), dim = p_AgD, dispstr = "un"),
    margins = c("norm","norm","binom",  "binom"),
    paramMargins = list(list(mean = m_AgD[1,1], sd = m_AgD[2,1]),
                        list(mean = m_AgD[1,2], sd = m_AgD[2,2]),
                        list(size = 1, prob = m_AgD[1,3]),
                        list(size = 1, prob = m_AgD[1,4]))
  )
mycov_AgD <-rMvdc(N_AgD, myMvd_AgD)

#Simulate binary outcome 
bX_AgD <- b[1]+b[-1]%*%t(mycov_AgD)
y_AgD <- rbinom(N_AgD,1,1/(1+exp(-bX_AgD)))

#Simulated IPD for the AgD study
AgD <- data.frame(cbind(y_AgD,mycov_AgD))
colnames(AgD) <- c("y",covariates_AgD)

#Generate summary statistics for the simulated AgD study
summary_AgD <- matrix(nrow=2,ncol=2+p_AgD)
colnames(summary_AgD) <-
  c("number of events", "number of patients",covariates_AgD)
rownames(summary_AgD)<-c("mean","sd")
summary_AgD[1,1:2] <-c(sum(AgD$y),N_AgD) 
summary_AgD[1,-(1:2)] <- colMeans(mycov_AgD)
summary_AgD[2,3:4] <- apply(mycov_AgD[,1:2],2,sd)

################################################################################################
#STEP5: Save data to be used in the STC
################################################################################################
#Data used in the STC related to IPD population
write.csv(IPD, "IPD_sim.csv", row.names = FALSE)

#Data used in the STC related to AgD population
write.csv(summary_AgD, "AgD_sim.csv")
