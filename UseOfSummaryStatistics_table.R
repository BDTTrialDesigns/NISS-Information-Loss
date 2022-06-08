# source("dataGenerationST2.R")
# setwd
#source("loadData.R")

# install.packages("devtools")
# library(devtools)
# install_github("starima74/AddInf", force=TRUE)
library(AddInf)

## here I redefine MMSE function to add bootstrap CI 
MMSE <- function (dd, theta.f, Add.Inf, nboots = 500, eig.cutoff = 0.9,  nsims = 1000, ...)  {
  library(MASS)
  p <- nrow(Add.Inf)
  n <- nrow(dd)
  thetahat <- theta.f(dd, ...)
  thetahat.length <- length(thetahat)
  betahat <- c()
  for (i in 1:p) {
    betahat.add <- Add.Inf$Functions[[i]](dd)
    betahat <- c(betahat, betahat.add)
  }
  betahat.length <- length(betahat)
  bootres <- matrix(NA, nrow = nboots, ncol = thetahat.length + 
                      betahat.length)
  for (i in 1:nboots) {
    ind <- sample(1:n, replace = T)
    boot.hat <- theta.f(dd[ind, ], ...)
    for (j in 1:p) boot.hat <- c(boot.hat, Add.Inf$Functions[[j]](dd[ind, 
    ]))
    bootres[i, ] <- boot.hat
  }
  K = var(bootres)
  ind1 <- 1:thetahat.length
  ind2 <- (thetahat.length + 1):(betahat.length + thetahat.length)
  K11 <- K[ind1, ind1]
  K12 <- K[ind1, ind2]
  K21 <- K[ind2, ind1]
  if (is.null(dim(bootres[, ind2]))) 
    Biases2 <- (mean(bootres[, ind2]) - unlist(Add.Inf$Means))^2
  if (!is.null(dim(bootres[, ind2]))) 
    Biases2 <- (colMeans(bootres[, ind2]) - unlist(Add.Inf$Means))^2
  Biases2 <- n * Biases2 %o% Biases2
  Biases2[Add.Inf$Biases == 0, ] <- 0
  Biases2[, Add.Inf$Biases == 0] <- 0
  K22 <- K[ind2, ind2] + Biases2
  K.add <- as.matrix(Matrix::bdiag(Add.Inf$Vars))
  betatilde <- unlist(Add.Inf$Means)
  SVD_var <- svd(K22 + K.add)
  ind.svd <- rep(1, length(SVD_var$d))
  if (length(SVD_var$d) > 1) 
    for (i in 2:length(SVD_var$d)) ind.svd[i] <- sum(SVD_var$d[1:(i - 
                                                                    1)])/sum(SVD_var$d) < eig.cutoff
  diag_elem <- rep(0, length(SVD_var$d))
  diag_elem[SVD_var$d > 0 & ind.svd == 1] <- 1/SVD_var$d[SVD_var$d > 
                                                           0 & ind.svd == 1]
  if (length(SVD_var$d > 0 & ind.svd == 1) == 1) 
    Kinv <- SVD_var$v %*% diag_elem %*% t(SVD_var$u)
  if (length(SVD_var$d > 0 & ind.svd == 1) > 1) 
    Kinv <- SVD_var$v %*% diag(diag_elem) %*% t(SVD_var$u)
  thetahat.MVAR = thetahat - K12 %*% Kinv %*% (betahat - betatilde)
  thetahat.MVAR.Var = K11 - K12 %*% Kinv %*% K21
  ### bootstrap
  thetahat.sims <- mvrnorm(nsims, mu = c(thetahat, betahat), 
                           Sigma = K)
  betatilde.sims <- mvrnorm(nsims, mu = betatilde, Sigma = K.add)
  thetaest.sims <- matrix(NA, nrow = nsims, ncol = length(thetahat))
  for (s in 1:nsims) {
    Biases2 <- (thetahat.sims[s, ind2] - betatilde.sims[s, 
    ])^2
    Biases2 <- n * Biases2 %o% Biases2
    Biases2[Add.Inf$Biases == 0, ] <- 0
    Biases2[, Add.Inf$Biases == 0] <- 0
    SVD_var <- svd(K[ind2, ind2] + Biases2 + K.add)
    ind.svd <- rep(1, length(SVD_var$d))
    if (length(SVD_var$d) > 1) 
      for (i in 2:length(SVD_var$d)) ind.svd[i] <- 
      sum(SVD_var$d[1:(i - 1)])/sum(SVD_var$d) < eig.cutoff
    diag_elem <- rep(0, length(SVD_var$d))
    diag_elem[SVD_var$d > 0 & ind.svd == 1] <- 1/SVD_var$d[SVD_var$d > 
                                                             0 & ind.svd == 1]
    if (length(SVD_var$d > 0 & ind.svd == 1) == 1) 
      Kinv <- SVD_var$v %*% diag_elem %*% t(SVD_var$u)
    if (length(SVD_var$d > 0 & ind.svd == 1) > 1) 
      Kinv <- SVD_var$v %*% diag(diag_elem) %*% t(SVD_var$u)
    thetaest.sims[s, ] = thetahat.sims[s, ind1] - K12 %*% 
      Kinv %*% (thetahat.sims[s, ind2] - betatilde.sims[s, 
      ])
  }
  ### bootstrap under theta = 0
  thetahat.sims0 <- mvrnorm(nsims, mu = c(0, betahat), 
                           Sigma = K)
  betatilde.sims0 <- mvrnorm(nsims, mu = betatilde, Sigma = K.add)
  thetaest.sims0 <- matrix(NA, nrow = nsims, ncol = length(thetahat))
  for (s in 1:nsims) {
    Biases2 <- (thetahat.sims0[s, ind2] - betatilde.sims0[s, 
    ])^2
    Biases2 <- n * Biases2 %o% Biases2
    Biases2[Add.Inf$Biases == 0, ] <- 0
    Biases2[, Add.Inf$Biases == 0] <- 0
    SVD_var <- svd(K[ind2, ind2] + Biases2 + K.add)
    ind.svd <- rep(1, length(SVD_var$d))
    if (length(SVD_var$d) > 1) 
      for (i in 2:length(SVD_var$d)) ind.svd[i] <- 
      sum(SVD_var$d[1:(i - 1)])/sum(SVD_var$d) < eig.cutoff
    diag_elem <- rep(0, length(SVD_var$d))
    diag_elem[SVD_var$d > 0 & ind.svd == 1] <- 1/SVD_var$d[SVD_var$d > 
                                                             0 & ind.svd == 1]
    if (length(SVD_var$d > 0 & ind.svd == 1) == 1) 
      Kinv <- SVD_var$v %*% diag_elem %*% t(SVD_var$u)
    if (length(SVD_var$d > 0 & ind.svd == 1) > 1) 
      Kinv <- SVD_var$v %*% diag(diag_elem) %*% t(SVD_var$u)
    thetaest.sims0[s, ] = thetahat.sims0[s, ind1] - K12 %*% 
      Kinv %*% (thetahat.sims0[s, ind2] - betatilde.sims0[s, 
      ])
  }
  list(Theta.Est = thetahat.MVAR, Theta.Hat = thetahat, Theta.Hat.Var = K11, 
       Theta.Est.sims = thetaest.sims, Theta.Est.sims0 = thetaest.sims0)
}

### end MMSE

load('TrialData.RData')

str(Sim$Sim_Datasets)

BMI0 <- Sim$Sim_Datasets$Main  #### main dataset
#BMI0 <- Sim$Sim_Datasets$Main2  #### main2 dataset

LM_z <- function(d) summary(lm(qnorm(BMI3) ~ R + covariate + qnorm(BMI1), data = d))$coef["R", c(1,2,4)]
LM_p <- function(d) summary(lm(BMI3 ~ R + covariate + BMI1, data = d))$coef["R", c(1,2,4)]

LM.int_z <- function(d) summary(lm(qnorm(BMI2) ~ R + covariate + qnorm(BMI1), data = d))$coef["R", c(1,2,4)]
LM.int_p <- function(d) summary(lm(BMI2 ~ R + covariate + BMI1, data = d))$coef["R", c(1,2,4)]

round(fit.p <- LM_p(BMI0),4)
round(c(fit.p[1]-1.96*fit.p[2],fit.p[1]+1.96*fit.p[2]),4)
round(fit.z <- LM_z(BMI0),4)
round(c(fit.z[1]-1.96*fit.z[2],fit.z[1]+1.96*fit.z[2]),4)

#############################################################3

theta.f.p <- function(d)  LM_p(d)[1]
theta.f.z <- function(d)  LM_z(d)[1]

Add.Info.Means <- list()
Add.Info.Vars <- list()
Add.Info.Functions <- list()
Add.Info.Biases <- list()

results <- data.frame(Method = c(rep("MVAR", 5),rep("MMSE",5)),
                      ExtData = c(1:5,1:5),
                      EffSE.p = rep(NA,10),
                      P.p = rep(NA,10),
                      CI.p = rep(NA,10),
                      EffSE.z = rep(NA,10),
                      P.z = rep(NA,10),
                      CI.z = rep(NA,10)
)

MMSE.CI.width <-1:10
for(i in 1:5) {
 if(i == 1) BMI1 <- Sim$Sim_Datasets$`Ext1-full`   # external with same sample size
 if(i == 2) BMI1 <- Sim$Sim_Datasets$`Ext2-double` # external with double sample size
 if(i == 3) BMI1 <- Sim$Sim_Datasets$`Ext3-half`   # external with half sample size
 if(i == 4) BMI1 <- Sim$Sim_Datasets$`Ext4-gdist`  # external with a strong departure from the null
 if(i == 5) BMI1 <- Sim$Sim_Datasets$`Ext5-mdist`  # external with a moderate departure from the null
 
 ### p chunk
 External.Info <- LM_p(BMI1)
 Add.Info.Means[[1]] <- External.Info[1]
 Add.Info.Vars[[1]] <- External.Info[2]^2
 Add.Info.Functions[[1]] <- function(d) LM_p(d)[1]

 Add.Info <- data.frame(Means = rep(NA, 1),
                       Vars = rep(NA, 1),
                       Functions = rep(NA, 1),
                       Biases = rep(NA, 1))
 Add.Info$Means = Add.Info.Means
 Add.Info$Vars = Add.Info.Vars
 Add.Info$Functions = Add.Info.Functions
 Add.Info$Biases = rep(NA, 1)

 set.seed(1234)
 res.MVAR <- MVAR(BMI0, theta.f.p, Add.Info, nboots = 500, eig.cutoff = 1)
 results[results$Method == "MVAR" & results$ExtData == i,3] <- 
   paste(round(res.MVAR$Theta.Est,4),"(", round(sqrt(res.MVAR$Theta.Est.Var),4),")",sep="")
 results[results$Method == "MVAR" & results$ExtData == i,4] <- 
   format(round(1-pchisq(res.MVAR$Theta.Est^2/res.MVAR$Theta.Est.Var,df=1),4), scientific=F)
 results[results$Method == "MVAR" & results$ExtData == i,5] <- 
   paste("[", round(res.MVAR$Theta.Est - 1.96 * sqrt(res.MVAR$Theta.Est.Var),4), ", ",
 round(res.MVAR$Theta.Est + 1.96 * sqrt(res.MVAR$Theta.Est.Var),4),"]",sep="")

 ### z chunk
 External.Info <- LM_z(BMI1)
 Add.Info.Means[[1]] <- External.Info[1]
 Add.Info.Vars[[1]] <- External.Info[2]^2
 Add.Info.Functions[[1]] <- function(d) LM_z(d)[1]
 
 Add.Info <- data.frame(Means = rep(NA, 1),
                        Vars = rep(NA, 1),
                        Functions = rep(NA, 1),
                        Biases = rep(NA, 1))
 Add.Info$Means = Add.Info.Means
 Add.Info$Vars = Add.Info.Vars
 Add.Info$Functions = Add.Info.Functions
 Add.Info$Biases = rep(NA, 1)
 
 set.seed(1234)
 res.MVAR <- MVAR(BMI0, theta.f.z, Add.Info, nboots = 500, eig.cutoff = 1)
 results[results$Method == "MVAR" & results$ExtData == i,6] <- 
   paste(round(res.MVAR$Theta.Est,4),"(", round(sqrt(res.MVAR$Theta.Est.Var),4),")",sep="")
 results[results$Method == "MVAR" & results$ExtData == i,7] <- 
   format(round(1-pchisq(res.MVAR$Theta.Est^2/res.MVAR$Theta.Est.Var,df=1),4), scientific=F)
 results[results$Method == "MVAR" & results$ExtData == i,8] <- 
   paste("[", round(res.MVAR$Theta.Est - 1.96 * sqrt(res.MVAR$Theta.Est.Var),4), ", ",
         round(res.MVAR$Theta.Est + 1.96 * sqrt(res.MVAR$Theta.Est.Var),4),"]",sep="")
 
 ### p chunk (MMSE)
 External.Info <- LM_p(BMI1)
 Add.Info.Means[[1]] <- External.Info[1]
 Add.Info.Vars[[1]] <- External.Info[2]^2
 Add.Info.Functions[[1]] <- function(d) LM_p(d)[1]
 
 Add.Info <- data.frame(Means = rep(NA, 1),
                        Vars = rep(NA, 1),
                        Functions = rep(NA, 1),
                        Biases = rep(NA, 1))
 Add.Info$Means = Add.Info.Means
 Add.Info$Vars = Add.Info.Vars
 Add.Info$Functions = Add.Info.Functions
 Add.Info$Biases = rep(1, 1)
 
 set.seed(1234)
 res.MMSE <- MMSE(BMI0, theta.f.p, Add.Info, nboots = 500, nsims = 10000, eig.cutoff = 1)
 results[results$Method == "MMSE" & results$ExtData == i,3] <- 
   paste(round(res.MMSE$Theta.Est,4),"(", round(sqrt(VAR <- var(res.MMSE$Theta.Est.sims)),4),")",sep="")
 results[results$Method == "MMSE" & results$ExtData == i,4] <- 
   format(round(mean(c(res.MMSE$Theta.Est.sims0)^2 > c(res.MMSE$Theta.Est)^2),4), scientific=F)
 results[results$Method == "MMSE" & results$ExtData == i,5] <- 
   paste("[", round(quantile(res.MMSE$Theta.Est.sims,prob=c(0.025)),4), ", ",
         round(quantile(res.MMSE$Theta.Est.sims,prob=c(0.975)),4),"]",sep="")
 
 ### z chank (MMSE)
 External.Info <- LM_z(BMI1)
 Add.Info.Means[[1]] <- External.Info[1]
 Add.Info.Vars[[1]] <- External.Info[2]^2
 Add.Info.Functions[[1]] <- function(d) LM_z(d)[1]
 
 Add.Info <- data.frame(Means = rep(NA, 1),
                        Vars = rep(NA, 1),
                        Functions = rep(NA, 1),
                        Biases = rep(NA, 1))
 Add.Info$Means = Add.Info.Means
 Add.Info$Vars = Add.Info.Vars
 Add.Info$Functions = Add.Info.Functions
 Add.Info$Biases = rep(1, 1)
 
 set.seed(1234)
 res.MMSE <- MMSE(BMI0, theta.f.z, Add.Info, nboots = 500, nsims = 10000, eig.cutoff = 1)
 results[results$Method == "MMSE" & results$ExtData == i,6] <- 
   paste(round(res.MMSE$Theta.Est,4),"(", round(sqrt(VAR <- var(res.MMSE$Theta.Est.sims)),4),")",sep="")
 results[results$Method == "MMSE" & results$ExtData == i,7] <- 
   format(round(mean(c(res.MMSE$Theta.Est.sims0)^2 > c(res.MMSE$Theta.Est)^2),4), scientific=F)
 results[results$Method == "MMSE" & results$ExtData == i,8] <- 
   paste("[", round(l.ci <- quantile(res.MMSE$Theta.Est.sims,prob=c(0.025)),4), ", ",
         round(u.ci <-quantile(res.MMSE$Theta.Est.sims,prob=c(0.975)),4),"]",sep="")
 MMSE.CI.width[i] <- u.ci - l.ci 
}

results

library(xtable)
xtable(results)
