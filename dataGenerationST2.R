#' Code from Sergey Tarima for generating BMI reduction datasets
#' This code generates the main dataset (BMI0) and 
#' 4 different external datasets (BMI1, BMI2, BM3 and BMI4)

library(MASS)
library(lme4)
set.seed(1234)

## GLOBAL VARIABLES

eff1 <- -0.015   # %BMI reduction at the early point (the intervention effect)
eff2 <- -0.0075  # %BMI reduction at the final endpoint (the intervention effect)
red1 <- -0.01     # %BMI reduction at the early point
red2 <- 0.00      # %BMI reduction at the final endpoint
lost1 <- 0.2   # probability of dropping prior to the early %BMI assessment
lost2 <- 0.4   # probability of dropping prior to the final %BMI assessment
baselineCovBMI1 <- -0.01 # the higher %BMI is the stronger %BMI reduction 
                         #is anticipated (early endpoint) toward %BMI = 0.5
baselineCovBMI2 <- -0.01 # the higher BMI is the stronger BMI reduction 
                         # is anticipated (final endpoint) toward %BMI = 0.5
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
SDs <- sqrt(c(0.02, 0.02, 0.02)) # the SDs are all 2% now
SIGMA <- diag(SDs) %*% CORR %*% diag(SDs)

############ FUNCTIONS (start)  ##############
############ FUNCTIONS (start)  ##############
############ FUNCTIONS (start)  ##############

## (1) data generating function
dataGeneration <- function(n = n0, 
                           seed = seed0, 
                           prob = prob0, 
                           baselineBMI = baselineBMI0,
                           baselineCovBMI = c(baselineCovBMI1, baselineCovBMI2),
                           additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
                           treatment_effect) {
  R <- rbinom(n, size = 1, prob = prob)
  covariate <- rnorm(n, mean = 1, sd = 1) 
  # note, it is better to pass the mean and SD of the covariate 
  # through the arguments; I am thinking that the covariate is AGE (?)
  # AGE is not assiciated with baseline, but maybe associated with 
  # the follow-up assessments
  MU <- data.frame(BMI1 = baselineBMI, 
                   BMI2 = baselineBMI + red1 + covariate*additionalCovBMI[1] + 
                     R*treatment_effect[1], 
                   BMI3 = baselineBMI + red2 + covariate*additionalCovBMI[2] + 
                     R*treatment_effect[2])
  BMI <- MASS::mvrnorm(n, mu = c(0,0,0), Sigma = SIGMA) + MU
  # below we add the effect of baseline on follow-up 
  # [the higher initial BMI the higher improvement we anticipate] 
  BMI[,-1] <- BMI[,-1] + (BMI[,1] - 0.5) %o% baselineCovBMI
  BMI$R <- R
  BMI[runif(n) < lost1, c("BMI2", "BMI3")] <- NA
  BMI[runif(n) < lost2, c("BMI3")] <- NA
  return(BMI)
}

## (2) plotting BMI data
BMIplot <- function(BMI, main = "BMI data") {
  n <- nrow(BMI)
  boxplot(c(BMI$BMI1, BMI$BMI2, BMI$BMI3)~
            c(rep(0, n), rep(12, n), rep(24, n)),
          xlab="Month", ylab="BMI", main = main)
  for(i in 1:10) {
    tmp <- BMI[i, ]
    points(c(rep(1, 1), rep(2, 1), rep(3, 1)),
           c(tmp$BMI1, tmp$BMI2, tmp$BMI3), pch = 19)
    lines(c(rep(1, 1), rep(2, 1), rep(3,1)),
          c(tmp$BMI1, tmp$BMI2, tmp$BMI3), lty = 2)
  }
}

############ FUNCTIONS (end) ##############
############ FUNCTIONS (end)  ##############
############ FUNCTIONS (end)  ##############


## DATASET 0 (main dataset)
n0 <- 452           # total sample size of the main study
seed0 <- 1234       # seed value to generate the main dataset
prob0 <- 0.5        # probability to assign the new treatment
baselineBMI0 <- 0.9 # mean baseline %BMI in main study 
eff10 <- -0.015   # %BMI reduction at the early time point (the intervention effect)
eff20 <- -0.0075  # %BMI reduction at the final endpoint (the intervention effect)

## DATASET 0.1 ('sister' trial)
n0 <- 200           # total sample size of the main study
seed01 <- 1235       # seed value to generate the main dataset
prob0 <- 0.5        # probability to assign the new treatment
baselineBMI01 <- 0.9 # mean baseline BMI in main study 
eff10 <- -0.015   # %BMI reduction at the early time point (the intervention effect)
eff20 <- -0.0075  # %BMI reduction at the final endpoint (the intervention effect)

## DATASET 1 (external dataset with a large sample)
n1 <- 2000          # total sample size of the 1st external dataset
seed1 <- 1236       # seed value to generate the 1st external dataset
prob1 <- 0.5        # probability to assign the new treatment
baselineBMI1 <- 0.9 # mean baseline %BMI in 1st external study
eff11 <- -0.015   # %BMI reduction at the early time point (the intervention effect)
eff21 <- -0.0075  # %BMI reduction at the final endpoint (the intervention effect)

## DATASET 2 (external dataset with a similar sample size)
n2 <- 300           # total sample size of the 1st external dataset
seed2 <- 1237       # seed value to generate the 1st external dataset
prob2 <- 0.5        # probability to assign the new treatment
baselineBMI2 <- 0.9 # mean baseline %BMI in 1st external study
eff12 <- -0.015   # %BMI reduction at the early time point (the intervention effect)
eff22 <- -0.0075  # %BMI reduction at the final endpoint (the intervention effect)

## DATASET 3 (external dataset with a small sample size)
n3 <- 50            # total sample size of the 1st external dataset
seed3 <- 1238       # seed value to generate the 1st external dataset
prob3 <- 0.5        # probability to assign the new treatment
baselineBMI3 <- 0.9 # mean baseline BMI in 1st external study
eff13 <- -0.015   # %BMI reduction at the early time point (the intervention effect)
eff23 <- -0.0075  # %BMI reduction at the final endpoint (the intervention effect)

#external data-sets - different treatment effect & baseline (same sample size)

## DATASET 4 (moderate conflict - treatment effect approx. 0.5 SD larger)
n4 <- 200          # total sample size of the 1st external dataset
seed4 <- 1239       # seed value to generate the 1st external dataset
prob4 <- 0.5        # probability to assign the new treatment
baselineBMI4 <- 0.9 # mean baseline BMI in 1st external study
eff14 <- eff10 + 0.5*0.013 
eff24 <- eff20 + 0.5*0.020 

## DATASET 5 (moderate conflict - treatment effect approx. 2SD larger)
n5 <- 200          # total sample size of the 1st external dataset
seed5 <- 1239       # seed value to generate the 1st external dataset
prob5 <- 0.5        # probability to assign the new treatment
baselineBMI5 <- 0.9 # mean baseline BMI in 1st external study
eff15 <- eff10 + 2*0.013 
eff25 <- eff20 + 2*0.020 

## DATA GENERATION
## DATA GENERATION
## DATA GENERATION

BMI0 <- dataGeneration(n0, 
                       seed = seed0, 
                       prob = prob0, 
                       baselineBMI = baselineBMI0,
                       baselineCovBMI = c(baselineCovBMI1,baselineCovBMI2),
                       additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
                       treatment_effect=c(eff10,eff20))

BMI0.1 <- dataGeneration(n0, 
                       seed = seed01, 
                       prob = prob0, 
                       baselineBMI = baselineBMI0,
                       baselineCovBMI = c(baselineCovBMI1,baselineCovBMI2),
                       additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
                       treatment_effect=c(eff10,eff20))

BMI1 <- dataGeneration(n1, 
                       seed=seed1, 
                       prob = prob1, 
                       baselineBMI = baselineBMI1,
                       baselineCovBMI = c(baselineCovBMI1,baselineCovBMI2),
                       additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
                       treatment_effect=c(eff11,eff21))

BMI2 <- dataGeneration(n2, 
                       seed=seed2, 
                       prob = prob2, 
                       baselineBMI = baselineBMI2,
                       baselineCovBMI = c(baselineCovBMI1,baselineCovBMI2),
                       additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
                       treatment_effect=c(eff12,eff22))

BMI3 <- dataGeneration(n3, 
                       seed=seed3, 
                       prob = prob3, 
                       baselineBMI = baselineBMI3,
                       baselineCovBMI = c(baselineCovBMI1,baselineCovBMI2),
                       additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
                       treatment_effect=c(eff13,eff23))

BMI4 <- dataGeneration(n4, 
                       seed=seed4, 
                       prob = prob4, 
                       baselineBMI = baselineBMI4,
                       baselineCovBMI = c(baselineCovBMI1,baselineCovBMI2),
                       additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
                       treatment_effect=c(eff14,eff24))

BMI5 <- dataGeneration(n5, 
                       seed=seed5, 
                       prob = prob5, 
                       baselineBMI = baselineBMI5,
                       baselineCovBMI = c(baselineCovBMI1,baselineCovBMI2),
                       additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
                       treatment_effect=c(eff15,eff25))

## (NOT RUN)
## plotting the generated BMI data
# BMIplot(BMI0, main = "Main BMI data")
# BMIplot(BMI1, main = "External large BMI data")
# BMIplot(BMI2, main = "External moderate BMI data")
# BMIplot(BMI3, main = "External small BMI data")
# BMIplot(BMI4, main = "External large BMI data\nbaseline BMI is similar to the main data")
