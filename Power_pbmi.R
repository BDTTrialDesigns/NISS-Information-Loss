####################################################
#
# Script to run the modified PowerPrior approach
# for linear regression
#
####################################################
library(rstan)
#library(brms)

rm(list=ls())
sm_lr <- stan_model(file="linear_reg_stan.stan", model_name='lralpha',verbose=FALSE)

hellinger_norm = function(x){
  # x = mu1, sigma1, mu2, sigma2
  sqrt(1 - sqrt(2*x[2]*x[4] / (x[2]^2 + x[4]^2))*exp(- ((x[1]-x[3])^2)/(4*(x[2]^2 + x[4]^2))) )
}

load('TrialData.RData')
th0=0
al=0.05

##### firs step, we check the distances between the main dataset and the external ones
BMI_r=list()
BMI_res=list()
fit_firststage = list()
for (i in 1:(length(Sim$Sim_Datasets)-1))
{
  BMI = Sim$Sim_Datasets[[i]]
  
  #only focus on late endpoint - discard any missings
  BMI_r[[i]]=BMI[apply(BMI,1,function(x) sum(is.na(x))==0),] #no NAs
  
  ####################################################################################################################################
  # run the model for the commensurability
  if(dim(BMI_r[[1]])[1] / dim(BMI_r[[i]])[1] <= 1){ # if the sample size of the historical datasets is higher than the actual one
    # we can use the results stored for the actual one in the first position of the list
    datastan = list(N = dim(BMI_r[[i]])[1],
                    Y = BMI_r[[i]]$BMI3,
                    K = 4,
                    X = cbind(Intercept = rep(1, dim(BMI_r[[i]])[1]), 
                              R = BMI_r[[i]]$R,
                              covariate = BMI_r[[i]]$covariate,
                              BMI1 = BMI_r[[i]]$BMI1),
                    alpha = dim(BMI_r[[1]])[1] / dim(BMI_r[[i]])[1],
                    priorMS = cbind(rep(0,4), rep(10,4)))
    
    set.seed(10)
    fit <- sampling(sm_lr, data = datastan, iter = 8000, 
                    chains = 4, control = list(adapt_delta = 0.85), cores=1)
    
    BMI_res[[i]] <- list(results = summary(fit, pars = c("b", "b_Intercept"), probs = c(0.025, 0.975))$summary[, c("mean", "sd")])
    dimnames(BMI_res[[i]]$results)[[1]] = c("R", "covariate", "BMI1", "Intercept")
    BMI_res[[i]]$dist = apply(cbind(BMI_res[[i]]$results,BMI_res[[1]]$results), 1, hellinger_norm)
  } else { # otherwise we have to run again the regression for the actual trial
    datastan = list(N = dim(BMI_r[[1]])[1],
                    Y = BMI_r[[1]]$BMI3,
                    K = 4,
                    X = cbind(Intercept = rep(1, dim(BMI_r[[1]])[1]), 
                              R = BMI_r[[1]]$R,
                              covariate = BMI_r[[1]]$covariate,
                              BMI1 = BMI_r[[1]]$BMI1),
                    alpha = dim(BMI_r[[i]])[1] / dim(BMI_r[[1]])[1],
                    priorMS = cbind(rep(0,4), rep(10,4)))
    
    set.seed(10)
    fit <- sampling(sm_lr, data = datastan, iter = 8000, 
                    chains = 4, control = list(adapt_delta = 0.85), cores=1)
    
    res1 <- summary(fit, pars = c("b", "b_Intercept"), probs = c(0.025, 0.975))$summary[, c("mean", "sd")]
    dimnames(res1)[[1]] = c("R", "covariate", "BMI1", "Intercept")
    
    datastan = list(N = dim(BMI_r[[i]])[1],
                    Y = BMI_r[[i]]$BMI3,
                    K = 4,
                    X = cbind(Intercept = rep(1, dim(BMI_r[[i]])[1]), 
                              R = BMI_r[[i]]$R,
                              covariate = BMI_r[[i]]$covariate,
                              BMI1 = BMI_r[[i]]$BMI1),
                    alpha = 1,
                    priorMS = cbind(rep(0,4), rep(10,4)))
    
    set.seed(10)
    fit <- sampling(sm_lr, data = datastan, iter = 8000, 
                    chains = 4, control = list(adapt_delta = 0.85), cores=1)
    
    BMI_res[[i]] <- list(results = summary(fit, pars = c("b", "b_Intercept"), probs = c(0.025, 0.975))$summary[, c("mean", "sd")])
    dimnames(BMI_res[[i]]$results)[[1]] = c("R", "covariate", "BMI1", "Intercept")
    BMI_res[[i]]$dist = apply(cbind(res1,BMI_res[[1]]$results), 1, hellinger_norm)
  }
  ####################################################################################################################################
  # full results of the trial using Bayesian regression
  datastan = list(N = dim(BMI_r[[i]])[1],
                  Y = BMI_r[[i]]$BMI3,
                  K = 4,
                  X = cbind(Intercept = rep(1, dim(BMI_r[[i]])[1]), 
                            R = BMI_r[[i]]$R,
                            covariate = BMI_r[[i]]$covariate,
                            BMI1 = BMI_r[[i]]$BMI1),
                  alpha = 1,
                  priorMS = cbind(rep(0,4), rep(10,4)))
  
  set.seed(10)
  fit <- sampling(sm_lr, data = datastan, iter = 8000, 
                  chains = 4, control = list(adapt_delta = 0.85), cores=1)
  
  fit_firststage[[i]] = fit
  
  BMI_res[[i]]$results_full = summary(fit, pars = c("b", "b_Intercept"), probs = c(0.025, 0.975))$summary[, c("mean", "sd", "2.5%", "97.5%")]
  dimnames(BMI_res[[i]]$results_full)[[1]] = c("R", "covariate", "BMI1", "Intercept")
  BMI_res[[i]]$pH0=1-pnorm(th0,BMI_res[[i]]$results_full[1,1],BMI_res[[i]]$results_full[1,2])
}


##### update results using the idea of power prior approach
# 1) ESS based on unit info variance / check variance and sd changing
# 2) discounting variance of (1-BMI_res[[i]]$dist)^2
BMI_borrow = list()
BMI_borrow[[1]] = BMI_res[[1]]
fit_full = list()
fit_full[[1]] = fit_firststage[[1]]
for (i in 2:(length(Sim$Sim_Datasets)-1) )
{
  unitinfo = (1/(BMI_res[[i]]$results_full[,2]^2) )/length(BMI_r[[i]]$BMI3)
  sigmas_new = (unitinfo*250)^(-1)
  datastan = list(N = dim(BMI_r[[1]])[1],
                  Y = BMI_r[[1]]$BMI3,
                  K = 4,
                  X = cbind(Intercept = rep(1, dim(BMI_r[[1]])[1]), 
                            R = BMI_r[[1]]$R,
                            covariate = BMI_r[[1]]$covariate,
                            BMI1 = BMI_r[[1]]$BMI1),
                  alpha = 1,
                  priorMS = cbind(BMI_res[[i]]$results_full[,1], 
                  #                BMI_res[[i]]$results_full[,2]/(1-BMI_res[[i]]$dist)^2))
                                 c(sqrt(sigmas_new[1]/(1-BMI_res[[i]]$dist[1])^2), rep(10,3)))
                                  #sigmas_new/(1-BMI_res[[i]]$dist)^2))
              )
  set.seed(10)
  fit <- sampling(sm_lr, data = datastan, iter = 8000, 
                  chains = 4, control = list(adapt_delta = 0.85), cores=1)
  fit_full[[i]] = fit
  BMI_borrow[[i]] = list(results=summary(fit, pars = c("b", "b_Intercept"), probs = c(0.025, 0.975))$summary[, c("mean", "sd", "2.5%", "97.5%")])
  dimnames(BMI_borrow[[i]]$results)[[1]] = c("R", "covariate", "BMI1", "Intercept")
  BMI_borrow[[i]]$pH0=1-pnorm(th0,BMI_borrow[[i]]$results[1,1],BMI_borrow[[i]]$results[1,2])
  BMI_borrow[[i]]$dist=BMI_res[[i]]$dist
}


################## check the results
# showing the results 
res_final = c(BMI_borrow[[1]]$results_full[1,], round(BMI_borrow[[1]]$pH0,4), BMI_borrow[[1]]$dist[1], original=BMI_borrow[[1]]$results_full[1,], 
              original.pHO=round(BMI_borrow[[1]]$pH0,4) )
for (i in 2:(length(Sim$Sim_Datasets)-1)){
  res_final = rbind(res_final, c(BMI_borrow[[i]]$results[1,], round(BMI_borrow[[i]]$pH0,4) , BMI_borrow[[i]]$dist[1], 
                                 BMI_res[[i]]$results_full[1,], round(BMI_res[[i]]$pH0,4 )) )
}
dimnames(res_final)[[1]] = c("Main", "Ext1-full", 
                             "Ext2-double", "Ext3-half","Ext4-gdist","Ext5-mdist")
dimnames(res_final)[[2]][c(5,6)] = c("pHO", "distance")
res_final
save.image(file="res_Power_pBMI.RData")



library(rtf)
library(packHV)
file_name="res_Power_pBMI"
tablerep = cbind(name=dimnames(res_final)[[1]], round(res_final,4)) 
rtffile <- RTF(paste(file_name,".doc", sep=""), width=11, height=8.5)  # this can be an .rtf or a .doc
addParagraph(rtffile, "This is the output \n")
addTable(rtffile, tablerep , col.justify=c("L", rep("R", dim(tablerep)[2]-1)),
         header.col.justify=c("L",rep("R", dim(tablerep)[2]-1)))
done(rtffile)



##### tentative of plot
dataplot = NULL
for (i in 1: (length(Sim$Sim_Datasets)-1)){
  param = extract(fit_full[[i]], pars = c("b", "b_Intercept"))
  dataplot = rbind(dataplot,
                   cbind(value = c(c(param$b), c(param$b_Intercept)),
                   pars = rep(c("R", "covariate", "BMI1", "Intercept"), each=dim(param$b)[1]),
                   scenario = rep(i, dim(param$b)[1]*4))
              )
}

dataplot = data.frame(dataplot)
str(dataplot)
dataplot$value = as.numeric(as.character(dataplot$value))

dplot2 = ggplot(dataplot[which(dataplot$pars=="R"),], aes(value)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line.x = element_line(colour = "black"),
        plot.title = element_text(size=20, face="bold", hjust = 0.5),
        axis.title.x = element_text(size=14, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_blank(),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(size=12, hjust = 0.5),
        legend.position = "top",
        legend.text = element_text(size=12)) +
  geom_density(aes(color=scenario), size=1.2) #+
  #facet_grid(. ~  pars)
#scale_x_continuous("minute",breaks=c(-90,0,90,200,400,600), limits=c(-100, 650)) + scale_y_continuous("", limit =c(0.5,1.5)) +

dplot2
ggsave("R_2_v4.jpeg", width = 25, height = 20, units = "cm") 

