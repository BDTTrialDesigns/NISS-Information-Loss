#New modified version of Sergey's code for the data-generating process

library(MASS)

## GLOBAL VARIABLES

prob0 <- 0.5       # probability to assign the new treatment
red1 <- -1     # BMI reduction at an early time point
red2 <- 0      # BMI reduction at the primary endpoint
baselineCovBMI1 <- -0.1  # the higher BMI is the stronger BMI reduction 
#is anticipated (early endpoint) toward BMI = 25
baselineCovBMI2 <- -0.05 # the higher BMI is the stronger BMI reduction 
# is anticipated (primary endpoint) toward BMI = 25
additionalCovBMI1 = 0 # for now  0
additionalCovBMI2 = 0 # for now  0

## note the CORR is not a correlation matrix
## SDs are not standard deviations
## this is just a way to create a covariance matrix
rho <- 0.75    # a parameter of a covariance matrix
rho23 <- 0.95  # a parameter of a covariance matrix
CORR  <- diag(rep(3,1))
CORR[CORR == 0] <- rho
CORR[2,3] <- CORR[3,2] <- rho23
SDs <- sqrt(c(2, 2, 2)) # the SDs are all 2% now
SIGMA <- diag(SDs) %*% CORR %*% diag(SDs)

############ FUNCTIONS (start)  ##############

## (1) data generating function
dataGeneration <- function(n = n0, 
                           seed = seed0, 
                           prob = prob0,
                           lost = c(0,0), 
                           baselineBMI = baselineBMI0,
                           baselineCovBMI = c(baselineCovBMI1, baselineCovBMI2),
                           additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
                           cov_par = c(1,1),
                           treatment_effect)
  {
  set.seed(seed) 
  R <- rbinom(n, size = 1, prob = prob)
  covariate <- rnorm(n, mean = cov_par[1], sd = cov_par[2])
  # MU <- data.frame(BMI1 = baselineBMI, 
  #                  BMI2 = baselineBMI + red1 + covariate*additionalCovBMI[1] + 
  #                    R*treatment_effect[1], 
  #                  BMI3 = baselineBMI + red2 + covariate*additionalCovBMI[2] + 
  #                    R*treatment_effect[2])
  # BMI <- MASS::mvrnorm(n, mu = c(0,0,0), Sigma = SIGMA) + MU
  #   # below we add the effect of baseline on follow-up 
  # [the higher initial BMI the higher improvement we anticipate] 
  # BMI[,-1] <- BMI[,-1] + (BMI[,1] - 25) %o% baselineCovBMI
  
  ##proposal:
  BMI1 = rnorm(n, baselineBMI, SDs[1])
  BMI2 = rnorm(n, BMI1  + (BMI1 - 25) * baselineCovBMI1 + red1 + covariate*additionalCovBMI[1] + 
                 R*treatment_effect[1], SDs[2])
  BMI3 = rnorm(n, BMI1  + (BMI1 - 25) * baselineCovBMI2 + red2 + covariate*additionalCovBMI[2] + 
                 R*treatment_effect[2], SDs[3])
  
  BMI = data.frame(BMI1=BMI1,BMI2=BMI2,BMI3=BMI3)
  BMI$R <- R
  BMI$covariate <- covariate
  if(lost[1]>0)
  BMI[sample(n,lost[1]), c("BMI2", "BMI3")] <- NA
  if(lost[2]>0)
  BMI[sample(n,lost[2]), c("BMI3")] <- NA
  return(BMI)
}

## (2) plotting BMI data
BMIplot <- function(BMI, main = "BMI data") {
  n <- nrow(BMI)
  boxplot(c(BMI$BMI1, BMI$BMI2, BMI$BMI3)~
            c(rep(0, n), rep(12, n), rep(24, n)),
          xlab="Month", ylab="BMI", main = main)
  for(i in 1:20) {
    tmp <- BMI[i, ]
    points(c(rep(1, 1), rep(2, 1), rep(3, 1)),
           c(tmp$BMI1, tmp$BMI2, tmp$BMI3), pch = 19)
    lines(c(rep(1, 1), rep(2, 1), rep(3,1)),
          c(tmp$BMI1, tmp$BMI2, tmp$BMI3), lty = 2)
  }
}


############ FUNCTIONS (end) ##############

n0 <- 452             #sample size of main study
baselineBMI0 <- 35    # mean baseline BMI in main study 
eff10 <- -1.5         # BMI reduction at an early time point (the intervention effect)
eff20 <- -0.75        # BMI reduction at the primary endpoint (the intervention effect)
lost10 <- 108         # number of drop-outs prior to the early BMI assessment
lost20 <- 366         # number of drop-outs prior to the primary BMI assessment
baselineBMI_ext <- 30 # baseline of external studies

#Datasets set-ups

settings=data.frame ('seed'=integer(),'n'=integer(),'baseline'=numeric(),'eff12'=numeric(),'eff24'=numeric(),'lost12'=integer(),'lost24'=integer())

settings['Main',]=  c(1234, n0, baselineBMI0, eff10, eff20, lost10, lost20)             
settings['Sister',]=c(1235, n0, baselineBMI0, eff10, eff20, lost10, lost20) 
settings['Ext1',]=  c(1236, n0, baselineBMI_ext, eff10, eff20, 0, 0) 
settings['Ext2',]=  c(1237, round(n0*10), baselineBMI_ext, eff10, eff20, 0, 0) 
settings['Ext3',]=  c(1238, round(n0/10), baselineBMI_ext, eff10, eff20, 0, 0) 
settings['Ext4',]=  c(1239, n0, baselineBMI_ext, eff10 + 0.075 , eff20 + 0.15, 0, 0) #approx 0.5 SDs away
settings['Ext5',]=  c(1240, n0, baselineBMI_ext, eff10 + 2*0.15 , eff20 + 2*0.30, 0, 0)  #approx 2 SDs away


## DATA GENERATION

BMI_data <- lapply(1:nrow(settings), function(x) 
                       dataGeneration(
                       n=settings$n[x], 
                       seed = settings$seed[x], 
                       prob=prob0,
                       lost= c(settings$lost12[x],settings$lost24[x]), 
                       baselineBMI = settings$baseline[x],
                       baselineCovBMI = c(baselineCovBMI1,baselineCovBMI2),
                       additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
                       treatment_effect=c(settings$eff12[x],settings$eff24[x])))

names(BMI_data) = rownames(settings)



## (NOT RUN)
## plotting the generated BMI data
# lapply(1:length(BMI_data), function(x) BMIplot(BMI_data[[x]], main = names(BMI_data)[x]))

