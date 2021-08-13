#' Code from Sergey for generating BMI reduction data
#' This code can be updated to generate external data in a similar manner
library(MASS)

n <- 200
set.seed(1234)
R <- rbinom(n,size = 1, prob=0.5)
eff1 <- -1.5   ## BMI reduction at an early  time point (the intervention effect)
eff2 <- -0.75 ## BMI reduction at the primary endpoint (the intervention effect)
red1 <- -1   ## BMI reduction at an early time point
red2 <- 0   ## BMI reduction at the primary endpoint
lost1 <- 0.2## probability of dropping prior to the early BMI assessment
lost2 <- 0.4 ## probability of dropping prior to the primary BMI assessment
baselineBMI <- 35

MU <- data.frame(BMI1=baselineBMI, 
                 BMI2=baselineBMI+red1 + R*eff1, 
                 BMI3=baselineBMI+red2 + R*eff2)
CORR  <- diag(rep(3,1))
CORR[CORR==0] <- 0.9
SDs <- sqrt(c(2,3,4))
SIGMA <- diag(SDs) %*% CORR %*% diag(SDs) 
BMI <- MASS::mvrnorm(n,mu=c(0,0,0),Sigma = SIGMA) + MU
BMI$R <- R
BMI[runif(n) < lost1, c("BMI2","BMI3")] <- NA
BMI[runif(n) < lost2, c("BMI3")] <- NA
colMeans(!is.na(BMI))

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
