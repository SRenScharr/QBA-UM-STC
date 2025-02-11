#Performing sensitivity analysis for unmeasured confounding using the simulated data
#Assuming covariate X4 was not measured in the AgD study

#Load library
library(dplyr)
library(copula)
library(ggplot2)

#Set working directory
#path <- ""
#setwd(path)

#Load the STC function
source("STC_boot.R")


#STEP1: read data 
#Define measured and unmeasured covariates
measured_cov <- c("X1","X2","X3")
unmeasured_cov <- "X4"

#Load the IPD study
IPD_data <- read.csv("IPD_sim.csv") 
r_B <- table(IPD_data$y)[2] #number of events
n_B <- nrow(IPD_data)

#Load the AgD study
AgD_data0 <- read.csv("AgD_sim.csv") 
AgD_data <- AgD_data0[,!names(AgD_data0) %in% unmeasured_cov] #only read X1-X3
r_A <- AgD_data$number.of.events[1]
n_A <- AgD_data$number.of.patients[1]


#STEP2: deriving the treatment effects using different methods
#Method 2.1: naive comparison without adjustment
p_A <- r_A/n_A
p_B <- r_B/n_B
var_A <- 1/(n_A*p_A*(1-p_A))
var_B <- 1/(n_B*p_B*(1-p_B));
log_OR_na <- log(p_B/(1-p_B))-log(p_A/(1-p_A))
log_OR_l_na <-  log_OR_na - 1.96*sqrt(var_B+var_A)
log_OR_U_na <-  log_OR_na + 1.96*sqrt(var_B+var_A)

OR_no_adjust <- exp(c(log_OR_na,log_OR_l_na,log_OR_U_na))
OR_no_adjust

#Method 2.2: STC adjusting only the measured/observed covariates
OR_pa_adjust <- STC_boot(IPD_data=IPD_data,
                         AgD_data=AgD_data,
                         measured_cov=measured_cov,
                         con_ind=c(1,2),#indicate that the first two covariates are continuous
                         n_b=1000,
                         n_me=1000) 
OR_pa_adjust


#Method 2.3: STC adjusting the unmeasured covariates assuming a fixed value for the marginal mean
#As an example, the mean of unmeasured X4 in the AgD study is assumed to be the same as the mean of X4 from the IPD study
OR_fu_adjust <- STC_boot(IPD_data=IPD_data,
                         AgD_data=AgD_data,
                         measured_cov=measured_cov,
                         unmeasured_cov=unmeasured_cov,
                         con_ind=c(1,2),
                         un_bin_mean=colMeans(IPD_data[unmeasured_cov]),
                         n_b=1000,
                         n_me=1000)
OR_fu_adjust


#STEP3: present the results in a table
round.n <- 4
OR_tab_META <- rbind(c("no adjustment", round(OR_no_adjust,round.n)),
                     c("adjusting only the observed", round(OR_pa_adjust,round.n)),
                     c("adjusting the unmeasured assuming a fixed value",round(OR_fu_adjust,round.n)))

op <- list(X=measured_cov,U=unmeasured_cov,OR_table=OR_tab_META)
op





#################################################################################################
#Investigate the impact of different values for the marginal mean of X4 in the AgD study using plots
#For illustration purposes, the marginal mean of X4 in both the AgD and IPD study are also evaluated
p_AgD <- AgD_data0[unmeasured_cov][1,]#the truth that has not been reported
p_IPD <- colMeans(IPD_data[unmeasured_cov])

p <- sort(unique(c(seq(0.05, 0.95, by = 0.05), p_AgD, p_IPD)))

#Calculate the estimate treatment effect given different values for the marginal mean of X4 in the AgD study
output_STC <- vector(mode = "list", length(p))
output_STC_adj <- array(NA, length(p))
output_STC_adj_ll <- array(NA, length(p))
output_STC_adj_ul <- array(NA, length(p))

for (i in 1:length(p)) {
  output_STC[[i]] <- STC_boot(IPD_data=IPD_data,
                              AgD_data=AgD_data,
                              measured_cov=measured_cov,
                              unmeasured_cov=unmeasured_cov,
                              con_ind=c(1,2),
                              un_bin_mean=p[i],
                              n_b=1000,
                              n_me=1000) #Note that the number of bootstrap is fixed using 1000 rather than 10000 to speedup the computation
  output_STC_adj[i] <- output_STC[[i]][1]
  output_STC_adj_ll[i] <- output_STC[[i]][2]
  output_STC_adj_ul[i] <- output_STC[[i]][3]
}


#plot
mydata <-
  data.frame(
    STC_OR = output_STC_adj,
    delta = p,
    STC_OR_ll = output_STC_adj_ll,
    STC_OR_ul = output_STC_adj_ul
  )

p_AgD_i <- which(mydata$delta==p_AgD)[1]
p_IPD_i <- which(mydata$delta==p_IPD)[1]

gp <- ggplot(data = mydata)+
  annotate("rect",ymin=-Inf,ymax=1,xmin=-Inf,xmax=Inf,fill="lightcoral",alpha = 0.4)+
  annotate("rect",ymin=1,ymax=Inf,xmin=-Inf,xmax=Inf,fill="darkseagreen3",alpha = 0.5)+
  geom_ribbon(aes(ymin = STC_OR_ll, ymax = STC_OR_ul, x = delta),alpha = 0.5) +
  geom_line(mapping = aes(x = delta, y = STC_OR),linewidth = 1.05) +
  scale_linetype_manual(name = "", values = c(2, 4), 
                        guide = guide_legend(override.aes = list(color = c("blue", "red"))))+
  geom_point(data = data.frame(x=c(p_IPD,p_AgD),
                               y=c(mydata$STC_OR[p_IPD_i],mydata$STC_OR[p_AgD_i]),
                               z=c("% of U in the IPD study","true % of U in the AgD study")),
             aes(x = x,y = y,shape = factor(z)),size = 3) +
  scale_y_continuous(name = "Adjusted odds ratio", breaks = c(0.5, 1, 1.5, 2, 2.5,3,3.5,4)) +
  scale_x_continuous(name = "Unmeasured covariate (U) in the AgD study: X4", breaks =
                       seq(0, 1, 0.1)) +
  #to plot the naive results
  geom_hline(aes(yintercept=OR_no_adjust[1], linetype="OR naive"), color = "red")+
  #to plot the partial adjustment results
  geom_hline(aes(yintercept=OR_pa_adjust[1], linetype="OR adjusting for observed covariates"), color = "blue")+
  labs(shape="") +
  theme(legend.position="top")

gp

ggsave( filename=paste("results.png",sep=''),
        plot=gp,
        width= 9, height = 6,units = "in",
        dpi=300
)

