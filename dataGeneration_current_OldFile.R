#New modified version of Sergey's code for the data-generating process

library(pwr)
library(MASS)
library(blockrand)
library(truncnorm)

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

## GLOBAL VARIABLES

prob0 <- 0.5      # probability to assign the new treatment
red1 <- -0.2      # %BMI reduction at the early point
red2 <- 0.00      # %BMI reduction at the final endpoint
additionalCovBMI1 = 0.1 #effect of covariate on the early endpoint
additionalCovBMI2 = 0.1 #effect of covariate on the final endpoint
cov_par = c(1,1)  #mean and sd of putative additional normal covariate

## CORR is a correlation matrix
## SDs are standard deviations
rho23 <- 0.9  # correlation between BMI2 and BMI3
CORR  <- diag(rep(2,1))
CORR[CORR == 0] <- rho23
sigma1 = 1
sigmas <- c(0.3, 0.3) # SDs
SIGMA <- diag(sigmas) %*% CORR %*% diag(sigmas)

rho23_1 <- 0.6  # low correlation between BMI2 and BMI3
CORR_1  <- diag(rep(2,1))
CORR_1[CORR_1 == 0] <- rho23_1
SIGMA1 <- diag(sigmas) %*% CORR_1 %*% diag(sigmas)

############ FUNCTIONS (start)  ##############

## (1) data generating function
dataGeneration <- function(n = n0, 
                           seed = seed0, 
                           lost = c(0,0), 
                           baselineBMI = baselineBMI0,
                           baselineCovBMI = c(baselineCovBMI1, baselineCovBMI2),
                           additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
                           cov_par = c(1,1),
                           treatment_effect = c(d,d),
                           Sigma)  
{ 
  set.seed(seed) 
  R <- as.numeric(as.character(blockrand(n, num.levels = 2, 
            levels = c(0,1), block.sizes = 2)$treatment[1:n])) 
  covariate <- rnorm(n, mean = cov_par[1], sd = cov_par[2]) 
  BMI1 = rtruncnorm(n, a=qnorm(0.85), mean = 0, sd = sigma1) 
  MU <- data.frame(BMI2 = BMI1 + red1 + 
      covariate * additionalCovBMI[1] + R * treatment_effect[1], 
                   BMI3 = BMI1 + red2 + 
      covariate * additionalCovBMI[2] + R * treatment_effect[2])
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

############ FUNCTIONS (end) ##############
d = 0.32
n0 <- 452               # sample size of main study
baselineBMI0 <- 1.5     # mean baseline BMI in main study (on Z scake)
eff10 <- -d/4-d/4         # BMI reduction at an early time point (the intervention effect)
eff20 <- -d/4             # BMI reduction at the primary endpoint (the intervention effect)
lost10 <- 250/2         # number of drop-outs prior to the early BMI assessment
lost20 <- 250           # number of drop-outs prior to the primary BMI assessment
baselineBMI_ext <- 1.5  # baseline of external studies

#Datasets set-ups

settings=data.frame ('seed'=integer(),
                     'n'=integer(),
                     'baseline'=numeric(),
                     'eff12'=numeric(),
                     'eff24'=numeric(),
                     'lost12'=integer(),
                     'lost24'=integer())

settings['Main',]        = c(10381, n0,          baselineBMI0,    eff10, eff20, lost10, lost20)             
settings['Sister',]      = c(1338, n0,          baselineBMI0,    eff10, eff20, lost10, lost20) 
settings['Ext1-full',]   = c(271,  n0,          baselineBMI_ext, eff10, eff20, 0, 0) 
settings['Ext2-double',] = c(277,  round(n0*2), baselineBMI_ext, eff10, eff20, 0, 0) 
settings['Ext3-half',]   = c(745,  round(n0/2), baselineBMI_ext, eff10, eff20, 0, 0) 

eff10 <- -d/4-d/4 - 0.2         # BMI reduction at an early time point (the intervention effect)
eff20 <- -d/4     - 0.2  
settings['Ext4-gdist',]=  c(274, n0, baselineBMI_ext, eff10, eff20, 0, 0) 

eff10 <- -d/4-d/4  -0.06        # BMI reduction at an early time point (the intervention effect)
eff20 <- -d/4      -0.06  
settings['Ext5-mdist',]=  c(389, n0, baselineBMI_ext, eff10, eff20, 0, 0)  

#add global variables

settings$'Pr_assignment_to_treatm'= rep(prob0, nrow(settings))
settings$'Percent_BMI_red_12'=rep(red1, nrow(settings))
settings$'Percent_BMI_red_24'=rep(red2, nrow(settings))
settings$'Covariate_mean'= rep(cov_par[1], nrow(settings))
settings$'Covariate_sd'= rep(cov_par[2], nrow(settings))
settings$'Effect_covariate_12'=rep(additionalCovBMI1, nrow(settings))
settings$'Effect_covariate_24'=rep(additionalCovBMI2, nrow(settings))
settings$'SD_baseline'=rep(sigmas[1], nrow(settings))
# settings$'SD_12'=rep(sigmas[2], nrow(settings))
# settings$'SD_24'=rep(sigmas[3], nrow(settings))
#settings$'rho_baseline.12'=rep(rho, nrow(settings))
#settings$'rho_baseline.24'=rep(rho, nrow(settings))
#settings$'rho_12.24'=rep(rho23, nrow(settings))

## DATA GENERATION

BMI_data <- lapply(1:nrow(settings), function(x) 
                       dataGeneration(
                       n=settings$n[x], 
                       seed = settings$seed[x], 
                       lost= c(settings$lost12[x],settings$lost24[x]), 
                       baselineBMI = settings$baseline[x],
                       baselineCovBMI = c(baselineCovBMI1,baselineCovBMI2),
                       additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
                       cov_par = cov_par,
                       treatment_effect=c(settings$eff12[x],settings$eff24[x]),
                       Sigma=SIGMA))
names(BMI_data) = rownames(settings) 

#### adding a dataset with small correlation

settings2=data.frame ('seed'=integer(),
                      'n'=integer(),
                      'baseline'=numeric(),
                      'eff12'=numeric(),
                      'eff24'=numeric(),
                      'lost12'=integer(),
                      'lost24'=integer())
eff10 <- -d/4-d/4         # BMI reduction at an early time point (the intervention effect)
eff20 <- -d/4             # BMI reduction at the primary endpoint (the intervention effect)
lost10 <- 250/2           # number of drop-outs prior to the early BMI assessment
lost20 <- 250  
rho23 <- 0.6  # correlation between BMI2 and BMI3
CORR  <- diag(rep(2,1))
CORR[CORR == 0] <- rho23
sigma1 = 1
sigmas <- c(0.3, 0.3) 
SIGMA1 <- diag(sigmas) %*% CORR %*% diag(sigmas)

#settings2['Main2',]=  c(5791, n0, baselineBMI0, eff10, eff20, lost10, lost20)   
settings2['Main2',]=  c(5792, n0, baselineBMI0, eff10, eff20, lost10, lost20)   
#add global variables

settings2$'Pr_assignment_to_treatm'= rep(prob0, nrow(settings2))
settings2$'Percent_BMI_red_12'=rep(red1, nrow(settings2))
settings2$'Percent_BMI_red_24'=rep(red2, nrow(settings2))
settings2$'Covariate_mean'= rep(cov_par[1], nrow(settings2))
settings2$'Covariate_sd'= rep(cov_par[2], nrow(settings2))
settings2$'Effect_covariate_12'=rep(additionalCovBMI1, nrow(settings2))
settings2$'Effect_covariate_24'=rep(additionalCovBMI2, nrow(settings2))
settings2$'SD_baseline'=rep(sigmas[1], nrow(settings2))

BMI_data2 <- lapply(1:nrow(settings2), function(x) 
  dataGeneration(
    n=settings2$n[x], 
    seed = settings2$seed[x], 
    lost= c(settings2$lost12[x],settings2$lost24[x]), 
    baselineBMI = settings2$baseline[x],
    baselineCovBMI = c(baselineCovBMI1,baselineCovBMI2),
    additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
    cov_par = cov_par,
    treatment_effect=c(settings2$eff12[x],settings2$eff24[x]),
    Sigma=SIGMA1))
names(BMI_data2) = rownames(settings2) 

Sim=list(Sim_Datasets=c(BMI_data, BMI_data2))
summary(lm(BMI3 ~ R + covariate + BMI1, data = Sim$Sim_Datasets$Main))
BMIplot(Sim$Sim_Datasets$Main[Sim$Sim_Datasets$Main$R == 1,])
BMIplot(Sim$Sim_Datasets$Main[Sim$Sim_Datasets$Main$R == 0,])
#setwd("/home/starima/Dropbox/tmp/Moreno")
save(Sim,file='TrialData_revisedST5.RData')

#load(file='TrialData_revisedST4.RData')

# checking datasets
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

library(xtable)
xtable(res_dataset, sanitize.text.function = identity, include.rownames=TRUE, digits=4)

# checking datasets
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




