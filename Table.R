#p scale
rm(list=ls())
load("~/Nextcloud/NISSworkingGroup/res_Power_pBMI.RData")
rm(list=setdiff(ls(), "res_final"))
load("~/Nextcloud/NISSworkingGroup/hierarchical_analysis_corr_mu_high_p.RData")
rm(list=setdiff(ls(), c("results_onestep","results_twostep","res_final")))

Power = round(res_final[,c(1,2,5,3,4)],4)
Hierarchical = round(t(results_onestep),4)
MAC =  round(t(results_twostep),4)

Bayes.results.p=data.frame('Approach'=c("Complete.data.only",rep(c("Hierarchical","Power","MAC"),each=5)),
                           'External.dataset'=c(rownames(Power)[1],rownames(Hierarchical)[1:5],rownames(Power)[2:6],rownames(MAC)[1:5]),
                           'Posterior.mean (SD)'=paste0(format(c(Power[1,1],Hierarchical[1:5,1],Power[2:6,1],MAC[1:5,1]),digits=2,scientific=FALSE)," (",format(c(Power[1,2],Hierarchical[1:5,2],Power[2:6,2],MAC[1:5,2]),digits=2,scientific=FALSE),")"),
                           'Posterior.pH0'=c(Power[1,3],Hierarchical[1:5,3],Power[2:6,3],MAC[1:5,3]),
                           '95CI'=paste0("[",format(c(Power[1,4],Hierarchical[1:5,4],Power[2:6,4],MAC[1:5,4]),digits=2,scientific=FALSE),",",format(c(Power[1,5],Hierarchical[1:5,5],Power[2:6,5],MAC[1:5,5]),digits=2,scientific=FALSE),"]"))


#z scale
rm(list=setdiff(ls(), "Bayes.results.p"))
load("~/Nextcloud/NISSworkingGroup/res_Power_zBMI.RData")
rm(list=setdiff(ls(), c("Bayes.results.p","res_final")))
load("~/Nextcloud/NISSworkingGroup/hierarchical_analysis_corr_mu_high_z.RData")
rm(list=setdiff(ls(), c("Bayes.results.p","results_onestep","results_twostep","res_final")))

Power = round(res_final[,c(1,2,5,3,4)],4)
Hierarchical = round(t(results_onestep),4)
MAC =  round(t(results_twostep),4)

Bayes.results.z=data.frame('Approach'=c("Complete.data.only",rep(c("Hierarchical","Power","MAC"),each=5)),
                           'External.dataset'=c(rownames(Power)[1],rownames(Hierarchical)[1:5],rownames(Power)[2:6],rownames(MAC)[1:5]),
                           'Posterior.mean (SD)'=paste0(format(c(Power[1,1],Hierarchical[1:5,1],Power[2:6,1],MAC[1:5,1]),digits=2,scientific=FALSE)," (",format(c(Power[1,2],Hierarchical[1:5,2],Power[2:6,2],MAC[1:5,2]),digits=2,scientific=FALSE),")"),
                           'Posterior.pH0'=c(Power[1,3],Hierarchical[1:5,3],Power[2:6,3],MAC[1:5,3]),
                           '95CI'=paste0("[",format(c(Power[1,4],Hierarchical[1:5,4],Power[2:6,4],MAC[1:5,4]),digits=2,scientific=FALSE),",",format(c(Power[1,5],Hierarchical[1:5,5],Power[2:6,5],MAC[1:5,5]),digits=2,scientific=FALSE),"]"))

Bayes.results=cbind.data.frame(Bayes.results.p,Bayes.results.z[,3:5])
Bayes.results$External.dataset=c('Main',rep(c('Full(NC)','Double(NC)','Half(NC)','Full(SC)','Full(MC)'),3))

library(xtable)
xtable(Bayes.results, digits = 4)