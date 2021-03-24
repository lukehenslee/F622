#==================================================================================================
# FISH 622 homework #1: Basic population dynamics and model fitting
# 2/3/2021
# Luke Henslee, College of Fisheries and Ocean Sciences, UAF
#
#==================================================================================================
#NOTES:

#==================================================================================================

# Load packages
library(tidyverse)
library(dplyer)
library(ggplot2)
library(bbmle)

# Set working directory
setwd("C:/Users/lhhenslee/Desktop/Luke/School/F622")
hw1.data <- file.path(getwd(), "data")
hw1.data


# Problem 2 ####

# 2.1- Load data set
troptuna <- read.csv(file=file.path(hw1.data, "troptuna.csv"))

# 2.2- Create prediction function
cpue_pred <- function(N0 = 10, r = 0.2, delta = 0.5, time = c(0:43)) {
  cpue <- N0*((1-delta)*exp(-1*r*time) + delta)
  return(cpue)
}


# 2.3- Create SSQ function
cpue_ssq <- function(N0=10, r=0.2, delta=0.1, time=c(0:43), obs.cpue=NULL) {
  pred.cpue <- cpue_pred(N0, r, delta, time)
  sq <- (obs.cpue - pred.cpue)^2
  ssq <- sum(sq)
  return(ssq)
}


# 2.4- Fit the model to the data
?mle2

obs.cpue <- troptuna$cpue

cpue_fit <- mle2(cpue_ssq, 
                 start=list(N0=12, r=0.1, delta=0.5),
                 data=list(time=c(0:43), obs.cpue=obs.cpue))

cpue_fit

coef(cpue_fit)


# 2.5- Create function to calculate true log-likelihood
?dnorm

cpue_normLike <- function(N0=10, r=0.2, delta=0.1, time=c(0:43), ln_sigma=1, cpue.obs=NULL) {
  sigma <- exp(ln_sigma)
  cpue.pred <- cpue_pred(N0, r, delta, time)
  NLL <- -1*sum(dnorm(cpue.obs, cpue.pred, sigma, log = TRUE))
  return(NLL)
}


# 2.6- Fit new model
cpue_fit_ML <- mle2(cpue_normLike, 
                 start=list(N0=12, r=0.1, delta=0.5, ln_sigma=1),
                 data=list(time=c(0:43), cpue.obs=obs.cpue),
                 method= "Nelder-Mead")

cpue_fit_ML

coef(cpue_fit_ML)

# 2.7- Explore confint() 
confint(cpue_fit)

confint(cpue_fit_ML)


# Problem 3 ####

# 3.1- Remove first 10 CPUE observations
obs.cpue.minus.10 <- troptuna[c(11:44),2]

cpue_fit_ML_minus_10 <- mle2(cpue_normLike,
                             start = list(N0=3, r=0.1, delta=0.5, ln_sigma=1),
                             data= list(time=c(0:33), cpue.obs = obs.cpue.minus.10))
cpue_fit_ML_minus_10

coef(cpue_fit_ML_minus_10)

confint(cpue_fit_ML_minus_10)

#3.2- Remove last 20 CPUE observations 
obs.cpue.minus.20 <- troptuna[c(1:24),2]

cpue_fit_ML_minus_20 <- mle2(cpue_normLike,
                             start = list(N0=12, r=0.1, delta=0.5, ln_sigma=1),
                             data = list(time=c(0:23), cpue.obs = obs.cpue.minus.20))

cpue_fit_ML_minus_20

coef(cpue_fit_ML_minus_20)

confint(cpue_fit_ML_minus_20)
