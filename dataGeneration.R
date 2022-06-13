#New modified version of Sergey's code for the data-generating process
library(pwr)
library(MASS)
library(blockrand)
library(truncnorm)
library(xtable)

############ FUNCTIONS  ##############

## (1) data generating function
dataGeneration <- function(treatment_effect,
                           rho23, 
                           n,
                           lost = c(0,0), 
                           seed, 
                           baselineBMI,
                           additionalCovBMI,
                           cov_par,
                           red,
                           sigmas)  
{ 
  set.seed(seed) 
  R <- as.numeric(as.character(blockrand(n, num.levels = 2, 
            levels = c(0,1), block.sizes = 2)$treatment[1:n])) 
  covariate <- rnorm(n, mean = cov_par[1], sd = cov_par[2]) 
  BMI1 = rtruncnorm(n, a=qnorm(0.85), mean = 0, sd = sigmas[1]) 
  MU <- data.frame(BMI2 = BMI1 + red[1] + 
      covariate * additionalCovBMI[1] + R * treatment_effect[1], 
                   BMI3 = BMI1 + red[2] + 
      covariate * additionalCovBMI[2] + R * treatment_effect[2])
  
  CORR  <- diag(rep(2,1)) # correlation matrix
  CORR[CORR == 0] <- rho23
  Sigma <- diag(sigmas[2:3]) %*% CORR %*% diag(sigmas[2:3])
  
  BMI <- MASS::mvrnorm(n, mu = c(0,0), Sigma = Sigma) + MU
  BMI$R <- R
  BMI$BMI1 <- pnorm(BMI1)
  BMI$BMI2 <- pnorm(BMI$BMI2)
  BMI$BMI3 <- pnorm(BMI$BMI3)
  BMI$covariate <- covariate
  if(lost[2]>0) BMI[(n-lost[2]+1):n, c("BMI3")] <- NA
  if(lost[1]>0) BMI[(n-lost[1]+1):n, c("BMI2", "BMI3")] <- NA
  return(BMI)
}

## (2) plotting BMI data
BMIplot <- function(BMI, main = "", ylab="",xlab="") {
  n <- nrow(BMI)
  boxplot(c(BMI$BMI1, BMI$BMI2, BMI$BMI3)~
            c(rep(0, n), rep(12, n), rep(24, n)),
          xlab="Month", main = main, range=0, ylab=ylab,
          ylim=c(0.5,1))
  for(i in ceiling(runif(20)*n)) {
    tmp <- BMI[i, ]
    points(c(rep(1, 1), rep(2, 1), rep(3, 1)),
           c(tmp$BMI1, tmp$BMI2, tmp$BMI3), pch = 19, cex =0.5)
    lines(c(rep(1, 1), rep(2, 1), rep(3,1)),
          c(tmp$BMI1, tmp$BMI2, tmp$BMI3), lty = 2)
  }
  points(c(rep(1, 1), rep(2, 1), rep(3, 1)),
         c(mean(BMI$BMI1,na.rm=TRUE), 
           mean(BMI$BMI2,na.rm=TRUE), 
           mean(BMI$BMI3,na.rm=TRUE)), pch = 19, cex = 2)
  lines(c(rep(1, 1), rep(2, 1), rep(3, 1)),
         c(mean(BMI$BMI1,na.rm=TRUE), 
           mean(BMI$BMI2,na.rm=TRUE), 
           mean(BMI$BMI3,na.rm=TRUE)), lwd = 3)
}


############ DATA GENERATION ##############

#Parmeters

settings=data.frame ('seed'=integer(),
                     'n'=integer(),
                     'eff12'=numeric(),
                     'eff24'=numeric(),
                     'lost12'=integer(),
                     'lost24'=integer(),
                     'rho23'=integer())

eff_size = pwr.t.test(n=452/2, 
                      power=0.9, 
                      type="two.sample", alternative = "two.sided")$d
d = round(eff_size, 2)

# Two-sample t test power calculation 
# 
# n = 226
# d = 0.3055905
# sig.level = 0.05
# power = 0.9
# alternative = two.sided
# 
# NOTE: n is number in *each* group

d = 0.32
n0 <- 452             # sample size of main study
eff10 <- -d/4-d/4     # BMI reduction at an early time point (the intervention effect)
eff20 <- -d/4         # BMI reduction at the primary endpoint (the intervention effect)
lost10 <- 250/2       # number of drop-outs prior to the early BMI assessment
lost20 <- 250         # number of drop-outs prior to the primary BMI assessment
rho23 <- 0.9         # correlation between the 12 and 24 month measurements (high correlation data-set)
rho23_1 <- 0.6         # correlation between the 12 and 24 month measurements (low correlation data-set)

settings['Main',]        = c(10381, n0,          eff10, eff20, lost10, lost20, rho23)             
settings['Ext1-full',]   = c(271,  n0,           eff10, eff20, 0, 0, rho23) 
settings['Ext2-double',] = c(277,  round(n0*2),  eff10, eff20, 0, 0, rho23) 
settings['Ext3-half',]   = c(745,  round(n0/2),  eff10, eff20, 0, 0, rho23) 

eff101 <- -d/4-d/4 - 0.2         # BMI reduction at an early time point (the intervention effect)
eff201 <- -d/4     - 0.2  
settings['Ext4-gdist',]=  c(274, n0,  eff101, eff201, 0, 0, rho23) 

eff102 <- -d/4-d/4  -0.06        # BMI reduction at an early time point (the intervention effect)
eff202 <- -d/4      -0.06  
settings['Ext5-mdist',]=  c(389, n0,  eff102, eff202, 0, 0, rho23)  

settings['Main2',]        = c(5792, n0,  eff10, eff20, lost10, lost20, rho23_1)             


#Global parameters

additionalCovBMI1 = 0.1 #effect of covariate on the early endpoint
additionalCovBMI2 = 0.1 #effect of covariate on the final endpoint
baselineBMI0 <- 1.5     # mean baseline BMI in main study (on Z scale)
red1 <- -0.2            # %BMI reduction at the early point
red2 <- 0.00            # %BMI reduction at the final endpoint
cov_par = c(1,1)        #mean and SD of putative additional normal covariate
sigmas = c(1, 0.3, 0.3) #SDs

settings$'BMI_red_12'=rep(red1, nrow(settings))
settings$'BMI_red_24'=rep(red2, nrow(settings))
settings$'baseline'=rep(baselineBMI0 , nrow(settings))
settings$'Covariate_mean'= rep(cov_par[1], nrow(settings))
settings$'Covariate_sd'= rep(cov_par[2], nrow(settings))
settings$'Effect_covariate_12'=rep(additionalCovBMI1, nrow(settings))
settings$'Effect_covariate_24'=rep(additionalCovBMI2, nrow(settings))
settings$'SD_baseline'=rep(sigmas[1], nrow(settings))
settings$'SD_12'=rep(sigmas[2], nrow(settings))
settings$'SD_24'=rep(sigmas[3], nrow(settings))

#Data generation

BMI_data <- lapply(1:nrow(settings), function(x) 
                       dataGeneration(
                         treatment_effect=c(settings$eff12[x],settings$eff24[x]),
                         rho23=settings$rho23[x],
                         n=settings$n[x],
                         lost = c(settings$lost12[x],settings$lost24[x]), 
                         seed = settings$seed[x],  
                         baselineBMI = settings$baseline[x],
                         additionalCovBMI=c(settings$Effect_covariate_12[x], settings$Effect_covariate_24[x]),
                         cov_par = c(settings$Covariate_mean[x], settings$Covariate_sd[x]),
                         red = c(settings$BMI_red_12[x], settings$BMI_red_24[x]),
                         sigmas = c(settings$SD_baseline[x], settings$SD_12[x], settings$SD_24[x])))
names(BMI_data) = rownames(settings) 


Sim=list(Sim_Datasets=BMI_data)

save(Sim,file='TrialData.RData')

############ PLOTS ##############

BMIplot(Sim$Sim_Datasets$Main[Sim$Sim_Datasets$Main$R == 1,])
BMIplot(Sim$Sim_Datasets$Main[Sim$Sim_Datasets$Main$R == 0,])

############ EFFECTS ##############

#Treatment effects on p scale
res_dataset = NULL
for (i in 1:length(Sim$Sim_Datasets)){
  BMI = Sim$Sim_Datasets[[i]]
  res2 = summary(lm(BMI3 ~ R + covariate + BMI1, data=BMI))$coefficients["R", c(1,2,4)]
  res_dataset = rbind(res_dataset, res2)
}
dimnames(res_dataset)[[1]] = names(Sim$Sim_Datasets) 
dimnames(res_dataset)[[2]] = c("Estimate", "Std. Error", "Pr(>|t|)") 
res_dataset[,3] = round(res_dataset[,3],4)
res_dataset

xtable(res_dataset, sanitize.text.function = identity, include.rownames=TRUE, digits=4)

#Treatment effects on z scale
res_dataset2 = NULL
for (i in 1:length(Sim$Sim_Datasets)){
  BMI = Sim$Sim_Datasets[[i]]
  res2 = summary(lm(qnorm(BMI3) ~ R + covariate + qnorm(BMI1), data=BMI))$coefficients["R", c(1, 2, 4)]
  res_dataset2 = rbind(res_dataset2, res2)
}
dimnames(res_dataset2)[[1]] = names(Sim$Sim_Datasets) 
dimnames(res_dataset2)[[2]] = c("Estimate", "Std. Error", "Pr(>|t|)") 
res_dataset2[,3] = round(res_dataset2[,3],4)
res_dataset2

xtable(res_dataset2, sanitize.text.function = identity, include.rownames=TRUE, digits=4)


############ FULL DATA (no missings) ##############

BMI_data_full <- lapply(c(1,7), function(x) 
  dataGeneration(
    treatment_effect=c(settings$eff12[x],settings$eff24[x]),
    rho23=settings$rho23[x],
    n=settings$n[x],
    lost = c(0,0), 
    seed = settings$seed[x],  
    baselineBMI = settings$baseline[x],
    additionalCovBMI=c(settings$Effect_covariate_12[x], settings$Effect_covariate_24[x]),
    cov_par = c(settings$Covariate_mean[x], settings$Covariate_sd[x]),
    red = c(settings$BMI_red_12[x], settings$BMI_red_24[x]),
    sigmas = c(settings$SD_baseline[x], settings$SD_12[x], settings$SD_24[x])))
names(BMI_data_full) = c('Main_full','Main2_full')

Sim_full=list(Sim_Datasets=BMI_data_full)

#Treatment effects on p scale
res_dataset_full = NULL
for (i in 1:length(Sim_full$Sim_Datasets)){
  BMI = Sim_full$Sim_Datasets[[i]]
  res2 = summary(lm(BMI3 ~ R + covariate + BMI1, data=BMI))$coefficients["R", c(1,2,4)]
  res_dataset_full = rbind(res_dataset_full, res2)
}
dimnames(res_dataset_full)[[1]] = names(Sim_full$Sim_Datasets) 
dimnames(res_dataset_full)[[2]] = c("Estimate", "Std. Error", "Pr(>|t|)") 
res_dataset_full[,3] = round(res_dataset_full[,3],4)
res_dataset_full

#Treatment effects on z scale
res_dataset_full2 = NULL
for (i in 1:length(Sim_full$Sim_Datasets)){
  BMI = Sim_full$Sim_Datasets[[i]]
  res2 = summary(lm(qnorm(BMI3) ~ R + covariate + qnorm(BMI1), data=BMI))$coefficients["R", c(1, 2, 4)]
  res_dataset_full2 = rbind(res_dataset_full2, res2)
}
dimnames(res_dataset_full2)[[1]] = names(Sim_full$Sim_Datasets) 
dimnames(res_dataset_full2)[[2]] = c("Estimate", "Std. Error", "Pr(>|t|)") 
res_dataset_full2[,3] = round(res_dataset_full2[,3],4)
res_dataset_full2



