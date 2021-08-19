#' Based on code from Sergey for generating BMI reduction data
#' Updated to generate external data in a similar manner
#' 
#' According to original study, children aged from 6 to 12; baseline reflects
#' BMI of children age 6 at 85% percentile.
library(MASS)

set.seed(1234)

###### Macro Variables ######
n <- 200
lower_age <- 6 
upper_age <- 12
eff_age <- 1 ## Slope for age effect on BMI
red1 <- -1   ## BMI reduction at an early time point
red2 <- 0   ## BMI reduction at the primary endpoint
lost1 <- 0.2## probability of dropping prior to the early BMI assessment
lost2 <- 0.4 ## probability of dropping prior to the primary BMI assessment
intercept_BMI <- 17

##### Functions #####
genData <- function(eff1, eff2) {
  R <- rbinom(n,size = 1, prob=0.5)
  age <- runif(n, lower_age, upper_age)
  baseline_BMI <- intercept_BMI + eff_age*age 
  
  MU <- data.frame(BMI1=baseline_BMI, 
                   BMI2=baseline_BMI + red1 + R*eff1, 
                   BMI3=baseline_BMI + red2 + R*eff2)
  CORR  <- diag(rep(3,1))
  CORR[CORR==0] <- 0.9
  SDs <- sqrt(c(2,3,4))
  SIGMA <- diag(SDs) %*% CORR %*% diag(SDs) 
  BMI <- MASS::mvrnorm(n,mu=c(0,0,0),Sigma = SIGMA) + MU
  BMI$R <- R
  BMI$age <- age
  BMI[runif(n) < lost1, c("BMI2","BMI3")] <- NA
  BMI[runif(n) < lost2, c("BMI3")] <- NA
  colMeans(!is.na(BMI))
  BMI
}

plotData <- function(BMI){
  ### plotting BMI data
  boxplot(c(BMI$BMI1,BMI$BMI2,BMI$BMI3)~
            c(rep(0,n),rep(12,n),rep(24,n)),
          xlab="Month", ylab="BMI", main = "BMI data")
  for(i in 1:10) {
    tmp <- BMI[i,]
    points(c(rep(1,1),rep(2,1),rep(3,1)),
           c(tmp$BMI1,tmp$BMI2,tmp$BMI3),pch=19)
    lines(c(rep(1,1),rep(2,1),rep(3,1)),
          c(tmp$BMI1,tmp$BMI2,tmp$BMI3),lty=2)
  }
}


## Main study + external (same)
eff1 <- -1.5 ## BMI reduction at an early  time point (the intervention effect)
eff2 <- -0.75 ## BMI reduction at the primary endpoint (the intervention effect)
BMI <- genData(eff1, eff2)
saveRDS(BMI, "data_main.rds")
plotData(BMI)

BMI <- genData(eff1, eff2)
saveRDS(BMI, "data_ext_same.rds")

## External study (mild difference; approx 1 std of treatment effect (0.16; 0.19)) 
BMI <- genData(eff1= -1.65, eff2= -0.95)
saveRDS(BMI, "data_ext_mild.rds")

## External study (large difference; approx 2 std)
BMI <- genData(eff1= -1.8, eff2= -1.15)
saveRDS(BMI, "data_ext_diff.rds")
