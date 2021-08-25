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
additionalCovBMI1 = 0 # for now  0
additionalCovBMI2 = 0 # for now  0

## CORR is a correlation matrix
## SDs are standard deviations
rho <- 0.75    # correlation between BMI1 and BMI3 (BMI1 and BMI2)
rho23 <- 0.95  # correlation between BMI2 and BMI3
CORR  <- diag(rep(3,1))
CORR[CORR == 0] <- rho
CORR[2,3] <- CORR[3,2] <- rho23
sigmas <- c(0.02, 0.02, 0.02) # all SDs are 2% now
SIGMA <- diag(sigmas) %*% CORR %*% diag(sigmas)

############ FUNCTIONS (start)  ##############
############ FUNCTIONS (start)  ##############
############ FUNCTIONS (start)  ##############

## (1) data generating function
dataGeneration <- function(n = n0, 
                           seed = seed0, 
                           prob = prob0, 
                           baselineBMI = baselineBMI0,
                           additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
                           treatment_effect, mz = NA, m = NA) {
  R <- rbinom(n, size = 1, prob = prob)
  covariate <- rnorm(n, mean = 1, sd = 1) 
  MU <- data.frame(BMI1 = baselineBMI, 
                   BMI2 = baselineBMI + red1 + 
                     covariate * additionalCovBMI[1] + R * treatment_effect[1], 
                   BMI3 = baselineBMI + red2 + 
                     covariate * additionalCovBMI[2] + R * treatment_effect[2])
  BMI <- MASS::mvrnorm(n, mu = c(0,0,0), Sigma = SIGMA) + MU
  BMI$R <- R
  if(is.na(m) | is.na(mz)) {
    BMI[runif(n) < lost1, c("BMI2", "BMI3")] <- NA
    BMI[runif(n) < lost2, c("BMI3")] <- NA
  }
  if(!is.na(m) & !is.na(mz)) {
    BMI[1:(n-mz), c("BMI2", "BMI3")] <- NA
    BMI[1:(n-m),  c("BMI3")] <- NA
  }
  return(BMI)
}

## (2) plotting BMI data
BMIplot <- function(BMI, main = "", ylab="Age and sex adjusted BMI percentile") {
  n <- nrow(BMI)
  boxplot(c(BMI$BMI1, BMI$BMI2, BMI$BMI3)~
            c(rep(0, n), rep(12, n), rep(24, n)),
          xlab="Month", main = main, range=0, ylab=ylab)
  for(i in ceiling(runif(10)*n)) {
    tmp <- BMI[i, ]
    points(c(rep(1, 1), rep(2, 1), rep(3, 1)),
           c(tmp$BMI1, tmp$BMI2, tmp$BMI3), pch = 19)
    lines(c(rep(1, 1), rep(2, 1), rep(3,1)),
          c(tmp$BMI1, tmp$BMI2, tmp$BMI3), lty = 2)
  }
}

############ FUNCTIONS (end) ##############
############ FUNCTIONS (end)  #############
############ FUNCTIONS (end)  #############


## DATASET 0 (main dataset)
n0 <- 452           # total sample size of the main study (pts with BMI1 observed)
seed0 <- 1234       # seed value to generate the main dataset
prob0 <- 0.5        # probability to assign the new treatment
baselineBMI0 <- 0.9 # mean baseline %BMI in main study 
eff10 <- -0.015   # %BMI reduction at the early time point (the intervention effect)
eff20 <- -0.0075  # %BMI reduction at the final endpoint (the intervention effect)
m <- 86 # number of patients with BMI1, BMI2 and BMI3 observed
mz <- 344 # number of patients with BMI1, BMI2 observed

## DATASET 0.1 ('sister' trial)
n01 <- 200           # total sample size of the main study
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
                       additionalCovBMI = c(additionalCovBMI1, additionalCovBMI2),
                       treatment_effect = c(eff10, eff20),
                       mz = mz, m = m)

BMI0.1 <- dataGeneration(n01, 
                       seed = seed01, 
                       prob = prob0, 
                       baselineBMI = baselineBMI0,
                       additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
                       treatment_effect=c(eff10,eff20))

BMI1 <- dataGeneration(n1, 
                       seed=seed1, 
                       prob = prob1, 
                       baselineBMI = baselineBMI1,
                       additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
                       treatment_effect=c(eff11,eff21))

BMI2 <- dataGeneration(n2, 
                       seed=seed2, 
                       prob = prob2, 
                       baselineBMI = baselineBMI2,
                       additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
                       treatment_effect=c(eff12,eff22))

BMI3 <- dataGeneration(n3, 
                       seed=seed3, 
                       prob = prob3, 
                       baselineBMI = baselineBMI3,
                       additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
                       treatment_effect=c(eff13,eff23))

BMI4 <- dataGeneration(n4, 
                       seed=seed4, 
                       prob = prob4, 
                       baselineBMI = baselineBMI4,
                       additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
                       treatment_effect=c(eff14,eff24))

BMI5 <- dataGeneration(n5, 
                       seed=seed5, 
                       prob = prob5, 
                       baselineBMI = baselineBMI5,
                       additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
                       treatment_effect=c(eff15,eff25))

## (NOT RUN)
## plotting the generated BMI data
# set.seed(123); BMIplot(BMI0, main = "Main BMI data")
# set.seed(123); BMIplot(BMI1, main = "External large BMI data")
# set.seed(123); BMIplot(BMI2, main = "External moderate BMI data")
# set.seed(123); BMIplot(BMI3, main = "External small BMI data")
