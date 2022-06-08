library(RBesT)
library(rstan)
#rstan_options(javascript=FALSE)
rm(list=ls())
out_dir='~/Nextcloud/NISSWorkingGroup/Files_for_repo'
setwd(out_dir)
corr= "high" #'low' # 
scale= "p" # "z" # 

load(paste0(out_dir,'/TrialData.RData'))

#High correlation main dataset
if (corr=='high')
  Sim$Sim_Datasets= Sim$Sim_Datasets[c(1,3:7)] 

#Low correlation main dataset
if (corr=='low')
  Sim$Sim_Datasets= Sim$Sim_Datasets[c(8,3:7)]

names(Sim$Sim_Datasets)=c('Main','Full(NC)','Double(NC)','Half(NC)','Full(SC)','Full(MC)')

################################################################################
## Summary statistics
################################################################################

th0=0
al=0.05

n_all=matrix(NA,length(Sim$Sim_Datasets),2)

summ_stats=rep(NA,length(Sim$Sim_Datasets))
sigmaj=rep(NA,length(Sim$Sim_Datasets))
sigma_unit=rep(NA,length(Sim$Sim_Datasets))
conf.int=matrix(NA,length(Sim$Sim_Datasets),2)

BMI_r=list()

for (i in 1:length(Sim$Sim_Datasets))
{
  BMI = Sim$Sim_Datasets[[i]]
  
  #only focus on late endpoint - discard any missings
  BMI_r[[i]]=BMI[apply(BMI,1,function(x) sum(is.na(x))==0),] #no NAs
  if(scale=="z"){
    BMI_r[[i]][,c('BMI1','BMI2','BMI3')]=qnorm(as.matrix(BMI_r[[i]][,c('BMI1','BMI2','BMI3')]))
  }
  n_all[i,]=c(sum(BMI_r[[i]]$R==0),sum(BMI_r[[i]]$R==1))
  
  lm.m=lm(BMI3 ~ R + covariate + BMI1, data=BMI_r[[i]])
  lmfit=summary(lm(BMI3 ~ R + covariate + BMI1, data=BMI_r[[i]]))
  summ_stats[i]= lmfit$coefficients['R',1]
  conf.int[i,]=confint(lm.m,'R', level = 0.95)
  sigmaj[i]=lmfit$coefficients['R',2]^2
  
  sigma_unit[i]=lmfit$coefficients['R',2]^2*sum(n_all[i,])
}

# frequentist confidence intervals for current and external

rownames(conf.int) = names(Sim$Sim_Datasets)
colnames(conf.int) = c('low','up')
names(summ_stats) = names(Sim$Sim_Datasets)

out=list('CI'=conf.int,'est'=summ_stats)
save(out,file=paste0(out_dir,'/95CI.freq.',scale,'.Corr_',corr,'.RData'))

# models
sm_lr <- stan_model(file="linear_reg_stan.stan", model_name='lralpha',verbose=FALSE)
sm_hier1 <- stan_model(file="ancova_hier_onestep_v2.stan", model_name='hier1',verbose=FALSE)
sm_hier2 <- stan_model(file="ancova_hier_twostep.stan", model_name='hier2',verbose=FALSE)

#########################################################################################################
## No borrowing
#########################################################################################################

datastan = list(N = dim(BMI_r[[1]])[1],
                Y = BMI_r[[1]]$BMI3,
                K = 4,
                X = cbind(Intercept = rep(1, dim(BMI_r[[1]])[1]), 
                          R = BMI_r[[1]]$R,
                          covariate = BMI_r[[1]]$covariate,
                          BMI1 = BMI_r[[1]]$BMI1),
                alpha = 1,
                priorMS = cbind(rep(0,4), rep(10,4)))

fit <- sampling(sm_lr, data = datastan, iter=8000, 
                chains = 4, control = list(adapt_delta = 0.85), cores=1, seed=10)

eff2 <- extract(fit, 'b', permuted = TRUE)$b[,1]#current trial estimate posterior samples
no.borrowing.stats=c('post.mean' = mean(eff2),'post.sd'=sd(eff2),'post_pr_nll'=mean(eff2>0) ,
                     'post.CI.low'= quantile(eff2,c(0.025)),'post.CI.up'= quantile(eff2,c(0.975)))
no.borrowing.post = fit

#########################################################################################################
## MAC - two-step
#########################################################################################################

ancova_hier_twostep = lapply(2:length(Sim$Sim_Datasets), function(z)
{
    datastan = list(M = 2,  
                    mean_stat = c(summ_stats[1],summ_stats[z]),
                    sd_stat = c(sqrt(sigmaj[1]),sqrt(sigmaj[z])),
                    sigma_nu=sigma_unit[z]/4)
    fit <- sampling(sm_hier2, data = datastan, iter=8000, 
                    chains = 4, control = list(adapt_delta = 0.99), cores=4, seed=10)
    
    eff2 <- extract(fit, 'b_treat', permuted = TRUE)$b_treat[,1] #current trial estimate posterior samples
    
    out.stats=c('post.mean' = mean(eff2),'post.sd'=sd(eff2),'post_pr_nll'=mean(eff2>0) ,
                'post.CI.low'= quantile(eff2,c(0.025)),'post.CI.up'= quantile(eff2,c(0.975)))
    out.post = fit
    
  out=list('stats'=out.stats,'samples'=out.post)
  out
})



results_twostep=data.frame(do.call(cbind,lapply(ancova_hier_twostep, function(z) z$'stats')))
rownames(results_twostep)=c('Post.mean','Post.sd','Post.Pr.Null','95.CI.low','95.CI.up')
colnames(results_twostep)=names(Sim$Sim_Datasets)[2:length(Sim$Sim_Datasets)]
results_twostep$'No borrowing'= no.borrowing.stats

#########################################################################################################
## Hierarchical model - one-step
#########################################################################################################

ancova_hier_onestep = lapply(2:length(Sim$Sim_Datasets), function(z)
{
    datastan = list(M = 2,
                    N = rowSums(n_all)[c(1,z)],
                    Y = BMI_r[[1]]$BMI3,
                    K = 4,
                    X = cbind(Intercept = rep(1, dim(BMI_r[[1]])[1]), 
                              R = BMI_r[[1]]$R,
                              covariate = BMI_r[[1]]$covariate,
                              BMI1 = BMI_r[[1]]$BMI1),
                    priorMS = cbind(rep(0,4), rep(10,4)),
                    Y0 = BMI_r[[z]]$BMI3,
                    X0 = cbind(Intercept = rep(1, dim(BMI_r[[z]])[1]), 
                               R = BMI_r[[z]]$R,
                               covariate = BMI_r[[z]]$covariate,
                               BMI1 = BMI_r[[z]]$BMI1),
                    sigma_nu=sigma_unit[z]/4)
    
    fit <- sampling(sm_hier1, data = datastan, iter=8000, 
                    chains = 4, control = list(adapt_delta = 0.99), cores=4, seed=10)
    
    
    eff2 <- extract(fit, 'b_treat', permuted = TRUE)$b_treat[,1] #current trial estimate posterior samples
    
    out.stats=c('post.mean' = mean(eff2),'post.sd'=sd(eff2),'post_pr_nll'=mean(eff2>0) ,
                'post.CI.low'= quantile(eff2,c(0.025)),'post.CI.up'= quantile(eff2,c(0.975)))
    out.post = fit
    
    out=list('stats'=out.stats,'samples'=out.post)
    out
})

##############################################################################################
##Results
##############################################################################################


results_onestep=data.frame(do.call(cbind,lapply(ancova_hier_onestep, function(z) z$'stats')))
rownames(results_onestep)=c('Post.mean','Post.sd','Post.Pr.Null','95.CI.low','95.CI.up')
colnames(results_onestep)=names(Sim$Sim_Datasets)[2:length(Sim$Sim_Datasets)]
results_onestep$'No borrowing'= no.borrowing.stats

save.image(file=paste0(out_dir,"/hierarchical_analysis_corr_mu_",corr,'_',scale,".RData"))

