#' Code from Sergey Tarima for generating BMI reduction datasets
#' This code generates the main dataset (BMI0) and 
#' 4 different external datasets (BMI1, BMI2, BM3 and BMI4)

library(MASS)
library(lme4)

############ FUNCTIONS (start)  ##############
############ FUNCTIONS (start)  ##############
############ FUNCTIONS (start)  ##############

## (1) data generating function
dataGeneration <- function(n = n0, 
                           seed = seed0, 
                           prob = prob0, 
                           baselineBMI = baselineBMI0,
                           baselineCovBMI = c(baselineCovBMI1, baselineCovBMI2)) {
  R <- rbinom(n, size = 1, prob = prob)
  MU <- data.frame(BMI1 = baselineBMI, 
                   BMI2 = baselineBMI + red1 + R*eff1, 
                   BMI3 = baselineBMI + red2 + R*eff2)
  CORR  <- diag(rep(3,1))
  CORR[CORR == 0] <- 0.9
  SDs <- sqrt(c(2, 3, 4))
  SIGMA <- diag(SDs) %*% CORR %*% diag(SDs) 
  BMI <- MASS::mvrnorm(n, mu = c(0,0,0), Sigma = SIGMA) + MU
  BMI[,-1] <- BMI[,-1] + (BMI[,1] - 25) %o% baselineCovBMI
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

## (3) summary of BMI reduction
## (this function generates some summary statistics on a BMI data)
SummaryBMIReduction <- function(BMI) {
  REDUCTION <- BMI[,2:3] - BMI[,1]
  DIFF = colMeans(REDUCTION[BMI$R == 1,], na.rm = TRUE) - 
    colMeans(REDUCTION[BMI$R == 0,], na.rm = TRUE) 
  VAR =  var(REDUCTION[BMI$R == 1,], na.rm = TRUE, 
             use = "pairwise.complete.obs") + 
         var(REDUCTION[BMI$R == 0,], na.rm = TRUE, 
             use = "pairwise.complete.obs") 
  return(list(DIFF = DIFF, 
              VAR = VAR, 
              n1 = sum(!is.na(REDUCTION[, 1])),
              n2 = sum(!is.na(REDUCTION[, 2]))
  ))
}

## (4) linear model for the early endpoint
SummaryBMI_LM1 <- function(BMI) {
  BMI$REDUCTION <- BMI[, 2] - BMI[, 1]
  S <- summary(lm(REDUCTION ~ BMI1 + R, data = BMI))$coef
  list(eff1 = S[3, 1], SE = S[3, 2])
}

### (5) linear model for the primary endpoint
SummaryBMI_LM2 <- function(BMI) {
  BMI$REDUCTION <- BMI[, 3] - BMI[, 1]
  S <- summary(lm(REDUCTION ~ BMI1 + R, data = BMI))$coef
  list(eff1 = S[3, 1], SE = S[3, 2])
}

## (6) a linear mixed linear model for both endpoints
SummaryBMI_LMER <- function(BMI) {
  n <- nrow(BMI)
  BMI_long <- data.frame(TIME = c(rep(0, n), rep(1, n)),
                    BMI1 = c(BMI$BMI1, BMI$BMI1),
                    ID = c(1:n, 1:n),
                    R = c(BMI$R, BMI$R),
                    REDUCTION = c(BMI[,2] - BMI[,1], BMI[,3] - BMI[,1])
  )
  BMI_long <- BMI_long[!is.na(BMI_long$REDUCTION),]
  f <- formula(REDUCTION ~ BMI1 + TIME + BMI1:TIME + 
                 I(R==1 & TIME == 0) + I(R==1 & TIME == 1) + (1|ID))
  fit <- lmer(f, data = BMI_long)
  list(ESTIMATE = summary(fit)$coef[,1], VCOV = vcov(fit))
}

############ FUNCTIONS (end) ##############
############ FUNCTIONS (end)  ##############
############ FUNCTIONS (end)  ##############

## GLOBAL VARIABLES

eff1 <- -1.5   # BMI reduction at an early time point (the intervention effect)
eff2 <- -0.75  # BMI reduction at the primary endpoint (the intervention effect)
red1 <- -1     # BMI reduction at an early time point
red2 <- 0      # BMI reduction at the primary endpoint
lost1 <- 0.2   # probability of dropping prior to the early BMI assessment
lost2 <- 0.4   # probability of dropping prior to the primary BMI assessment
baselineCovBMI1 <- -0.1  # the higher BMI is the stronger BMI reduction 
                         #is anticipated (early endpoint) toward BMI = 25
baselineCovBMI2 <- -0.05 # the higher BMI is the stronger BMI reduction 
                         # is anticipated (primary endpoint) toward BMI = 25

## DATASET 0 (main dataset)
n0 <- 200          # total sample size of the main study
seed0 <- 1234      # seed value to generate the main dataset
prob0 <- 0.5       # probability to assign the new treatment
baselineBMI0 <- 35 # mean baseline BMI in main study 

## DATASET 1 (external dataset with a large sample)
n1 <- 2000         # total sample size of the 1st external dataset
seed1 <- 1235      # seed value to generate the 1st external dataset
prob1 <- 0.5       # probability to assign the new treatment
baselineBMI1 <- 30 # mean baseline BMI in 1st external study

## DATASET 2 (external dataset with a similar sample size)
n2 <- 300          # total sample size of the 1st external dataset
seed2 <- 1236      # seed value to generate the 1st external dataset
prob2 <- 0.5       # probability to assign the new treatment
baselineBMI2 <- 27 # mean baseline BMI in 1st external study

## DATASET 3 (external dataset with a small sample size)
n3 <- 50           # total sample size of the 1st external dataset
seed3 <- 1237      # seed value to generate the 1st external dataset
prob3 <- 0.5       # probability to assign the new treatment
baselineBMI3 <- 30 # mean baseline BMI in 1st external study

## DATASET 4 (external dataset with a large sample; similar baseline BMI)
n4 <- 2000         # total sample size of the 1st external dataset
seed4 <- 1238      # seed value to generate the 1st external dataset
prob4 <- 0.5       # probability to assign the new treatment
baselineBMI4 <- 35 # mean baseline BMI in 1st external study

## DATA GENERATION
## DATA GENERATION
## DATA GENERATION

BMI0 <- dataGeneration(n0, 
                       seed = seed0, 
                       prob = prob0, 
                       baselineBMI = baselineBMI0,
                       baselineCovBMI = c(baselineCovBMI1,baselineCovBMI2))

BMI1 <- dataGeneration(n1, 
                       seed=seed1, 
                       prob = prob1, 
                       baselineBMI = baselineBMI1,
                       baselineCovBMI = c(baselineCovBMI1,baselineCovBMI2))

BMI2 <- dataGeneration(n2, 
                       seed=seed2, 
                       prob = prob2, 
                       baselineBMI = baselineBMI2,
                       baselineCovBMI = c(baselineCovBMI1,baselineCovBMI2))

BMI3 <- dataGeneration(n3, 
                       seed=seed3, 
                       prob = prob3, 
                       baselineBMI = baselineBMI3,
                       baselineCovBMI = c(baselineCovBMI1,baselineCovBMI2))

BMI4 <- dataGeneration(n4, 
                       seed=seed4, 
                       prob = prob4, 
                       baselineBMI = baselineBMI4,
                       baselineCovBMI = c(baselineCovBMI1,baselineCovBMI2))

## (NOT RUN)
## plotting the generated BMI data
# BMIplot(BMI0, main = "Main BMI data")
# BMIplot(BMI1, main = "External large BMI data")
# BMIplot(BMI2, main = "External moderate BMI data")
# BMIplot(BMI3, main = "External small BMI data")
# BMIplot(BMI4, main = "External large BMI data\nbaseline BMI is similar to the main data")

## calculating some summary statistics [BMI reduction]
# SummaryBMIReduction(BMI0)
# SummaryBMIReduction(BMI1)
# SummaryBMIReduction(BMI2)
# SummaryBMIReduction(BMI3)
# SummaryBMIReduction(BMI4)

## fitting linear and a linear mixed models
# SummaryBMI_LM1(BMI1)
# SummaryBMI_LM2(BMI1)
# SummaryBMI_LMER(BMI1)
