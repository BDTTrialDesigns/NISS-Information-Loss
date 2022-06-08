source("dataGenerationST2.R")

####### FUNCTIONS (start) ########
####### FUNCTIONS (start) ########
####### FUNCTIONS (start) ########

## (1) summary of BMI reduction
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

## (2) linear model for the early endpoint
SummaryBMI_LM1 <- function(BMI) {
  BMI$REDUCTION <- BMI[, 2] - BMI[, 1]
  S <- summary(lm(REDUCTION ~ BMI1 + R, data = BMI))$coef
  list(eff = S[3, 1], SE = S[3, 2])
}

### (3) linear model for the primary endpoint
SummaryBMI_LM2 <- function(BMI) {
  BMI$REDUCTION <- BMI[, 3] - BMI[, 1]
  S <- summary(lm(REDUCTION ~ BMI1 + R, data = BMI))$coef
  list(eff = S[3, 1], SE = S[3, 2])
}

## (4) a linear mixed linear model for both endpoints
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

###### FUNCTIONS (END) #######
###### FUNCTIONS (END) #######
###### FUNCTIONS (END) #######

####################### EXAMPLE 1 ####################################
###### (Use of External Information on early endpoint) ###############
###### (External Information is available via a summary statistic) ###
####################### EXAMPLE 1 ####################################

## Suppose BMI1 was previously analyzed b a group of researchers.
## They have applied a ANCOVA model to predict 12-month change in BMI.
## In their model, authors modelled the effect of the intervention (R)
## controlling for baseline BMI. In their publication, they reported 
## the adjusted for baseline BMI effect of intervention R on 12-month BMI change.

SummaryBMI_LM1(BMI1)
#$eff
#[1] -0.02220584
#
#$SE
#[1] 0.004514242

## Note that in the original dataset BMI1 is not available and 24-month 
## intervention is not reported.

## Our objective is to estimate 24-month effect of the intervention R adjusted 
## for baseline BMI.

## We can reproduce the (adjusted) effect of R on 12-month BMI change.

SummaryBMI_LM1(BMI0)
# $eff
# [1] -0.03939954
#
# $SE
# [1] 0.01422625

## The 24-month effect can be estimated on our dataset using another ANCOVA

SummaryBMI_LM2(BMI0)
# $eff
# [1] -0.04097859
# 
# $SE
# [1] 0.02024492

## We will use "AddInf" package to incorporate this external information
## in estimation of 24-month change in BMI adjusted for baseline BMI

## (uncomment "devtools" installation if needed)
## install.packages("devtools")
library(devtools)
install_github("starima74/AddInf", force=TRUE)
library(AddInf)

## First, we define a function which calculates the 24-month adjusted 
## for baseline effect.

theta.f <- function(d) {
  SummaryBMI_LM2(d)$eff
}

## Define components of a structure where the additional information is stored.
Add.Info.Means <- list()
Add.Info.Vars <- list()
Add.Info.Functions <- list()
Add.Info.Biases <- list()

## fill in information from the external data source + function used to 
## calculate the (adjusted) effect at 12 months

External.Info <- SummaryBMI_LM1(BMI1)

Add.Info.Means[[1]] <- External.Info$eff
Add.Info.Vars[[1]] <- External.Info$SE^2
Add.Info.Functions[[1]] <- function(d) SummaryBMI_LM1(d)$eff

## Create a data frame with external information
Add.Info <- data.frame(Means = rep(NA, 1),
                       Vars = rep(NA, 1),
                       Functions = rep(NA, 1),
                       Biases = rep(NA, 1))
Add.Info$Means = Add.Info.Means
Add.Info$Vars = Add.Info.Vars
Add.Info$Functions = Add.Info.Functions
Add.Info$Biases = rep(NA, 1)


####################### EXAMPLE 1.1 ####################################
###### (apply minimum variance estimation with external information) ###
####################### EXAMPLE 1.1 ####################################

set.seed(1234)
res.MVAR <- MVAR(BMI0, theta.f, Add.Info, nboots = 500, eig.cutoff = 1)

str(res.MVAR)

# $ Theta.Est    : num [1, 1] -0.0254
# $ Theta.Est.Var: num [1, 1] 0.000234
# $ Theta.Hat    : num -0.041
# $ Theta.Hat.Var: num 0.000432

# asymptotic 95% confidence interval of the improved estimate
c(res.MVAR$Theta.Est - 1.96 * sqrt(res.MVAR$Theta.Est.Var), res.MVAR$Theta.Est + 1.96 * sqrt(res.MVAR$Theta.Est.Var))

#          2.5%        97.5% 
#[1] -0.055382892  0.004528637

# asymptotic 95% confidence interval of the empirical estimate based on main the data only
c(res.MVAR$Theta.Hat - 1.96 * sqrt(res.MVAR$Theta.Hat.Var), res.MVAR$Theta.Hat + 1.96 * sqrt(res.MVAR$Theta.Hat.Var))

#          2.5%        97.5% 
# [1] -0.0817165690 -0.0002406048

####################### EXAMPLE 1.2 ###############################
###### (apply minimum MSE estimation with external information) ###
####################### EXAMPLE 1.2 ###############################

# Now we assume that the external data may come from a slightly different 
# data source, or the intervention of this another dataset was slightly different
# Then, we apply minimum MSE approach which is robust to the departures from 
# the null model. The robustness is achieved by suppressing the effect of external
# data when the external data cannot be used to estimate the effect of interest

## Create a data frame with external information
Add.Info <- data.frame(Means = rep(NA, 1),
                       Vars = rep(NA, 1),
                       Functions = rep(NA, 1),
                       Biases = rep(NA, 1))
Add.Info$Means = Add.Info.Means
Add.Info$Vars = Add.Info.Vars
Add.Info$Functions = Add.Info.Functions
Add.Info$Biases = rep(1, 1) # here we say that the bias in estimation of 
                            # the effect of interest is possible

res.MMSE <- MMSE(BMI0, theta.f, Add.Info, nboots = 500, eig.cutoff = 1)
str(res.MMSE)

# $ Theta.Est     : num [1, 1] -0.0266
# $ Theta.Hat     : num -0.041
# $ Theta.Hat.Var : num 0.000406
# $ Theta.Est.sims: num [1:1000, 1:2000] -0.0228 -0.0254 -0.0434 -0.029 -0.058 ...

# asymptotic 95% confidence interval
quantile(res.MMSE$Theta.Est.sims,prob=c(0.025,0.975))

#          2.5%        97.5% 
# -0.0688934453 -0.0003629102 

# the MMSE function also allows to perform MVAR estimation

Add.Info$Biases = rep(0, 1) # here we say that the bias in estimation of 
                            # the effect of interest is NOT possible

res.MMSE2 <- MMSE(BMI0, theta.f, Add.Info, nboots = 500, eig.cutoff = 1)
str(res.MMSE2)

# $ Theta.Est     : num [1, 1] -0.0268
# $ Theta.Hat     : num -0.041
# $ Theta.Hat.Var : num 0.000374
# $ Theta.Est.sims: num [1:1000, 1:2000] -0.00461 -0.03756 -0.02416 -0.0413 -0.03456 ...

# asymptotic 95% confidence interval
quantile(res.MMSE2$Theta.Est.sims,prob=c(0.025,0.975))
#          2.5%        97.5% 
#  -0.055979938  0.003857322 