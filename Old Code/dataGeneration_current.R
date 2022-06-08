#New modified version of Sergey's code for the data-generating process

library(MASS)

## GLOBAL VARIABLES

prob0 <- 0.5       # probability to assign the new treatment
red1 <- -0.01     # %BMI reduction at the early point
red2 <- 0.00      # %BMI reduction at the final endpoint
additionalCovBMI1 = 0.005 #effect of being a monoparental household on the early endpoint
additionalCovBMI2 = 0.005 #effect of being a monoparental household on the final endpoint
cov_prob  = 0.2   #probability of monoparental household

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

## (1) data generating function
dataGeneration <- function(n = n0, 
                           seed = seed0, 
                           prob = prob0,
                           lost = c(0,0), 
                           baselineBMI = baselineBMI0,
                           baselineCovBMI = c(baselineCovBMI1, baselineCovBMI2),
                           additionalCovBMI=c(additionalCovBMI1, additionalCovBMI2),
                           cov_prob = cov_prob,
                           treatment_effect)
  {
  set.seed(seed) 
  R <- rbinom(n, size = 1, prob = prob)
  covariate <- rbinom(n, size = 1, prob = cov_prob) 
  MU <- data.frame(BMI1 = baselineBMI, 
                   BMI2 = baselineBMI + red1 + 
                     covariate * additionalCovBMI[1] + R * treatment_effect[1], 
                   BMI3 = baselineBMI + red2 + 
                     covariate * additionalCovBMI[2] + R * treatment_effect[2])
  BMI <- MASS::mvrnorm(n, mu = c(0,0,0), Sigma = SIGMA) + MU
  BMI$R <- R
  BMI$monoparental <- covariate
  if(lost[2]>0)
  BMI[(n-lost[2]+1):n, c("BMI3")] <- NA
  if(lost[1]>0)
  BMI[(n-lost[1]+1):n, c("BMI2", "BMI3")] <- NA
 
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

n0 <- 452               # sample size of main study
baselineBMI0 <- 0.9     # mean baseline BMI in main study 
eff10 <- -0.015         # BMI reduction at an early time point (the intervention effect)
eff20 <- -0.0075        # BMI reduction at the primary endpoint (the intervention effect)
lost10 <- 108           # number of drop-outs prior to the early BMI assessment
lost20 <- 366           # number of drop-outs prior to the primary BMI assessment
baselineBMI_ext <- 0.9  # baseline of external studies

#Datasets set-ups

settings=data.frame ('seed'=integer(),'n'=integer(),'baseline'=numeric(),'eff12'=numeric(),'eff24'=numeric(),'lost12'=integer(),'lost24'=integer())

settings['Main',]=  c(1234, n0, baselineBMI0, eff10, eff20, lost10, lost20)             
settings['Sister',]=c(1235, n0, baselineBMI0, eff10, eff20, lost10, lost20) 
settings['Ext1',]=  c(1236, n0, baselineBMI_ext, eff10, eff20, 0, 0) 
settings['Ext2',]=  c(1237, round(n0*100), baselineBMI_ext, eff10, eff20, 0, 0) 
settings['Ext3',]=  c(1238, round(n0/10), baselineBMI_ext, eff10, eff20, 0, 0) 
settings['Ext4',]=  c(1239, n0, baselineBMI_ext, eff10 + 0.5*0.0015 , eff20 + 0.5*0.003, 0, 0) #approx 0.5 SDs away
settings['Ext5',]=  c(1240, n0, baselineBMI_ext, eff10 + 3*0.0015 , eff20 + 3*0.003, 0, 0)  #approx 3 SDs away

#add global variables

settings$'Pr_assignment_to_treatm'= rep(prob0, nrow(settings))
settings$'Percent_BMI_red_12'=rep(red1, nrow(settings))
settings$'Percent_BMI_red_24'=rep(red2, nrow(settings))
settings$'Pr_monoparental'= rep(cov_prob, nrow(settings))
settings$'Effect_monoparental_12'=rep(additionalCovBMI1, nrow(settings))
settings$'Effect_monoparental_12'=rep(additionalCovBMI2, nrow(settings))
settings$'SD_baseline'=rep(sigmas[1], nrow(settings))
settings$'SD_12'=rep(sigmas[2], nrow(settings))
settings$'SD_24'=rep(sigmas[3], nrow(settings))
settings$'rho_baseline.12'=rep(rho, nrow(settings))
settings$'rho_baseline.24'=rep(rho, nrow(settings))
settings$'rho_12.24'=rep(rho23, nrow(settings))

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
                       cov_prob = cov_prob,
                       treatment_effect=c(settings$eff12[x],settings$eff24[x])))

names(BMI_data) = rownames(settings)

Sim=list(Sim_Datasets=BMI_data,Sim_Parameters=settings)
save(Sim,file='~/TrialData.RData')

## (NOT RUN)
## plotting the generated BMI data

# set.seed(12345)
# cairo_pdf(file = '~/Documents/NISSWorkingGroup/main.pdf', width=6, height=6, pointsize=9)
# BMIplot(BMI_data[[1]], main = names(BMI_data)[1])
# graphics.off()

# lapply(1:length(BMI_data), function(x) BMIplot(BMI_data[[x]], main = names(BMI_data)[x]))


