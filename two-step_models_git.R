library(RBesT)

######## Models for primary endpoint only

# IPD for external studies may or may not be available - either way, we assume that a
# summary statistics for the treatment effect and its SE is available (two-step approach)

rm(list=ls())

load('~/Documents/NISSWorkingGroup/TrialData.RData')

n_tot=matrix(NA,length(Sim$Sim_Datasets),2)
n_all=matrix(NA,length(Sim$Sim_Datasets),2)

rownames(n_all) = rownames(n_tot) = names(Sim$Sim_Datasets)
colnames(n_all) = colnames(n_tot) = c('Control', 'Treated')

summ_stats=rep(NA,length(Sim$Sim_Datasets))
sigmaj=rep(NA,length(Sim$Sim_Datasets))
BMI_r=list()

for (i in 1:length(Sim$Sim_Datasets))
{
  BMI = Sim$Sim_Datasets[[i]]
  
  n_tot[i,]=c(sum(BMI$R==0),sum(BMI$R==1))
  
  #only focus on late endpoint - discard any missings
  BMI_r[[i]]=BMI[apply(BMI,1,function(x) sum(is.na(x))==0),] #no NAs
  n_all[i,]=c(sum(BMI_r[[i]]$R==0),sum(BMI_r[[i]]$R==1))
  
  lmfit=summary(lm(BMI3 ~ R + covariate + BMI1, data=BMI_r[[i]]))
  summ_stats[i]= lmfit$coefficients['R',1]
  sigmaj[i]=lmfit$coefficients['R',2]^2
  
  if(i==1)
    sigma_unit=lmfit$coefficients['R',2]^2*sum(n_all[1,])
}

############################################################################################
##analysis

th0=0
al=0.05

summ_stats_borrowing_external = lapply(3:length(Sim$Sim_Datasets), function(z)
{
  #No borrowing
  
  post.var_NB=  sigmaj[1]
  post.mean_NB = summ_stats[1]
  post.CI_NB = c(post.mean_NB + qnorm(al/2) * sqrt(post.var_NB),post.mean_NB - qnorm(al/2) * sqrt(post.var_NB))
  post_pr_nll_NB=1-pnorm(th0,post.mean_NB,sqrt(post.var_NB))
  
  #Full borrowing
  
  post.var_FB=   1/(1/sigmaj[z]+1/sigmaj[1])
  post.mean_FB = post.var_FB*(summ_stats[z]/sigmaj[z] + summ_stats[1]/sigmaj[1])
  post.CI_FB = c(post.mean_FB + qnorm(al/2) * sqrt(post.var_FB),post.mean_FB - qnorm(al/2) * sqrt(post.var_FB))
  post_pr_nll_FB= 1-pnorm(th0,post.mean_FB,sqrt(post.var_FB))
  
  
  #EB power prior (for normal case, equivalent to EB commensurate)
  
  priorPars=c(summ_stats[z],sigmaj[z])
  d = priorPars[2] / (pmax((summ_stats[1] - priorPars[1]) ^ 2, sigmaj[1] + priorPars[2]) - sigmaj[1])
  
  prior.sd.EB = sqrt(priorPars[2]/d)
  post.var_EB=  1/(1/prior.sd.EB^2+1/sigmaj[1])
  post.mean_EB = post.var_EB*(priorPars[1]/prior.sd.EB^2 + summ_stats[1]/sigmaj[1])
  post.CI_EB = c(post.mean_EB + qnorm(al/2) * sqrt(post.var_EB),post.mean_EB - qnorm(al/2) * sqrt(post.var_EB))
  post_pr_nll_EB=1-pnorm(th0,post.mean_EB,sqrt(post.var_EB))
  
  #robust mixture prior (RBesT package)
  
  rob_prior <- mixnorm(c(0.8,summ_stats[z],sqrt(sigmaj[z])),c(0.2,summ_stats[z],sqrt(sigma_unit)), param = 'ms')
  post_mix <- postmix(rob_prior, m=summ_stats[1], se=sqrt(sigmaj[1]))
  
  post.mean_Mix = post_mix[1,]%*%post_mix[2,]
  post.var_Mix = post_mix[1,]%*% (post_mix[3,]^2  + post_mix[2,]^2) -post.mean_Mix^2 
  post_pr_nll_Mix <- pmix(post_mix, 0, lower.tail = FALSE)
  post.CI_Mix = c(qmix(post_mix, al/2),qmix(post_mix, 1-al/2))
  
  
  out=matrix(c(
    post.mean_NB,post.var_NB,post.CI_NB,post_pr_nll_NB,
    post.mean_EB,post.var_EB,post.CI_EB,post_pr_nll_EB,
    post.mean_Mix,post.var_Mix,post.CI_Mix,post_pr_nll_Mix,
    post.mean_FB,post.var_FB,post.CI_FB,post_pr_nll_FB),5,4)
  
  colnames(out)=c('No.borrowing','EB.Power/Comm','Robust.Mixture','Full.borrowing')
  rownames(out)=c('Post.mean','Post.var','95.CI.low','95.CI.up','Post.Pr.Null')
  
  out
  
})

names(summ_stats_borrowing_external)=names(Sim$Sim_Datasets)[3:length(Sim$Sim_Datasets)]

