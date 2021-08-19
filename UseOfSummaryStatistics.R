source("dataGenerationST.R")

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
# $eff
# [1] -1.525917
# 
# $SE
# [1] 0.0566189

## Note that in the original dataset BMI1 is not available and 24-month 
## intervention is not reported.

## Our objective is to estimate 24-month effect of the intervention R adjusted 
## for baseline BMI.

## We can reproduce the (adjusted) effect of R on 12-month BMI change.

SummaryBMI_LM1(BMI0)
# $eff
# [1] -1.670822
# 
# $SE
# [1] 0.190133

## The 24-month effect can be estimated on our dataset using another ANCOVA

SummaryBMI_LM2(BMI0)
# $eff
# [1] -1.084476
# 
# $SE
# [1] 0.2687321

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
                       Biases = rep(NA, ))
Add.Info$Means = Add.Info.Means
Add.Info$Vars = Add.Info.Vars
Add.Info$Functions = Add.Info.Functions

## apply minimum variance estimation with external information

set.seed(1234)
res <- MVAR(BMI0, theta.f, Add.Info, nboots = 500, eig.cutoff = 1)

res

# $Theta.Est
#            [,1]
# [1,] -0.9679861
#
# $Theta.Est.Var
#            [,1]
# [1,] 0.04400151
#
# $Theta.Hat
# [1] -1.084476
#
# $Theta.Hat.Var
# [1] 0.07194423

