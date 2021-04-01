#==================================================================================================
# FISH 622 Homework #3: Growth and per-recruit analysis
# 3/23/2021
# Luke Henslee, College of Fisheries and Ocean Sciences, UAF
#
#==================================================================================================
#NOTES:

#==================================================================================================

# Load packages and set working directory ####
library(tidyverse)
library(ggplot2)
library(bbmle)
library(dplyr)
library(manipulate)

# Set working directory
setwd("C:/Users/lhhenslee/Desktop/Luke/School/F622/data")

# Problem 1- Estimating differences in Northern Rockfish growth ####

# 1.1- Load dataset 
goa <- read.csv("goa_race_specimen.csv", header=TRUE, skip=6)
str(goa)
head(goa)

# 1.2- Filter data
nr.dat <- goa %>% filter(Common.Name=="northern rockfish",
                         !is.na(Age..years.), !is.na(Length..mm.),
                         Sex!="Unknown")

# 1.3- Create predictive function for length using LVB model
pred_LVB <- function(age, Linf, k, t0) {
  Lt <- Linf * (1 - exp(-1 * k * (age - t0)))
  return(Lt)
}

# 1.4- Create NLL function for LVB model
NLL_LVB <- function(ln_Linf, ln_k, t0, ln_sigma, obs.age, obs.length) {
  
  # Exponentiate parameters
  Linf <- exp(ln_Linf)
  k <- exp(ln_k)
  sigma <- exp(ln_sigma)
  
  # Calculate predictions under LVB model
  pred.length <- pred_LVB(obs.age, Linf, k, t0)
  
  # Calculate log likelihood using dnorm()
  NLL <- -1 * dnorm(x=log(obs.length+1e-6),
                    mean=log(pred.length+1e-6), sd=sigma, log=TRUE)
  
  # Calculate total negative log likelihood
  return(sum(NLL))
}

# 1.5 & 1.6- Fit model to data for female and male Northern Rockfish and extract
# parameter values

  # Females
nr.female <- filter(nr.dat, nr.dat$Sex == "Female")

LVB_fit_female <- mle2(NLL_LVB,
                       start = list(ln_Linf = log(max(nr.female$Length..mm.)), 
                                    ln_k = log(0.2), t0 = 0, ln_sigma = log(0.5)),
                       data = list(obs.age = nr.female$Age..years., 
                                   obs.length = nr.female$Length..mm.),
                       method = "Nelder-Mead",
                       optimizer = "nlminb",
                       control = list(maxit = 1e6))

Linf.female <- as.numeric(exp(coef(LVB_fit_female)[1]))
k.female <- as.numeric(exp(coef(LVB_fit_female)[2]))
t0.female <- as.numeric(coef(LVB_fit_female)[3])

  # Males
nr.male <- filter(nr.dat, nr.dat$Sex == "Male")

LVB_fit_male <- mle2(NLL_LVB,
                       start = list(ln_Linf = log(max(nr.male$Length..mm.)), 
                                    ln_k = log(0.2), t0 = 0, ln_sigma = log(0.5)),
                       data = list(obs.age = nr.male$Age..years., 
                                   obs.length = nr.male$Length..mm.),
                       method = "Nelder-Mead",
                       optimizer = "nlminb",
                       control = list(maxit = 1e6))

Linf.male <- as.numeric(exp(coef(LVB_fit_male)[1]))
k.male <- as.numeric(exp(coef(LVB_fit_male)[2]))
t0.male <- as.numeric(coef(LVB_fit_male)[3])

# 1.7- Report parameter estimates
  # For females
cat(paste("LVB parameter estimates for female northern rockfish\n", 
          "Linf estimate:", exp(coef(LVB_fit_female)[1]), "\n", "k estimate:", 
          exp(coef(LVB_fit_female)[2]), "\n", "t0 estimate:", coef(LVB_fit_female)[3]), 
    "\n", "sigma estimate:", exp(coef(LVB_fit_female)[4]))

  # For males
cat(paste("LVB parameter estimates for male northern rockfish\n", 
          "Linf estimate:", exp(coef(LVB_fit_male)[1]), "\n", "k estimate:", 
          exp(coef(LVB_fit_male)[2]), "\n", "t0 estimate:", coef(LVB_fit_male)[3]), 
    "\n", "sigma estimate:", exp(coef(LVB_fit_male)[4]))

summary(LVB_fit_female)
summary(LVB_fit_male)

# 1.8- Plot model output versus observations
  # Females
nr.female$LVB <- pred_LVB(nr.female$Age..years., Linf.female, k.female, t0.female)

ggplot(nr.female, aes(x = Age..years.)) +
  geom_line(aes(y = LVB), color = "red", size = 1.5) +
  geom_point(aes(y = Length..mm.)) + 
  ylim(0, 500) +
  labs(title = "Observed (points) and predicted (line) lengths at age for female
  Northern Rockfish", x = "Age", y = "Length (mm)")

  # Males  
nr.male$LVB <- pred_LVB(nr.male$Age..years., Linf.male, k.male, t0.male)

ggplot(nr.male, aes(x = Age..years.)) +
  geom_line(aes(y = LVB), color = "blue", size = 1.5) +
  geom_point(aes(y = Length..mm.)) + 
  ylim(0, 500) +
  labs(title = "Observed (points) and predicted (line) lengths at age for male
  Northern Rockfish", x = "Age", y = "Length (mm)")

# Problem 3- Pollock per-recruit analysis####

# 3.1- Load pollock data
pol.dat <- read.csv(file="pollock_race_specimen.csv", header=TRUE, skip=6)

  # Subset for complete observations
pol.dat <- pol.dat %>% filter(!is.na(Age..years.),
                              !is.na(Length..mm.),
                              !is.na(Weight..gm.))

# 3.2- Create predictive function for weight using allometric weight-length relationship
pred_wl <- function(length, alpha, beta) {
    wl <- alpha * (length ^ beta)
    return(wl)
}

#3.3- Create NLL function
NLL_wl <- function(ln_alpha, ln_beta, ln_sigma, obs.length, obs.weight) {
  # Exponentiate parameters
  alpha <- exp(ln_alpha)
  beta <- exp(ln_beta)
  sigma <- exp(ln_sigma)
  
  # Create model predictions
  pred.weight <- pred_wl(obs.length, alpha, beta)
  
  # Calculate log-likelihood (lognormal distribution)
  logLike <- dnorm(x=log(obs.weight), mean=log(pred.weight), sd=sigma, log=TRUE)
  
  # Calculate log likelihood using dnorm()
  NLL <- -1 * dnorm(x = log(obs.weight + 1e-6),
                    mean = log(pred.weight + 1e-6), sd = sigma, log=TRUE)
  
  # Calculate total negative log likelihood
  return(sum(NLL))
}

# 3.4- Fit model
wl_fit <- mle2(NLL_wl,
               start = list(ln_alpha=log(1e-6), ln_beta=log(3), 
                            ln_sigma=log(0.2)),
               data = list(obs.length = pol.dat$Length..mm., 
                           obs.weight = pol.dat$Weight..gm.),
               method = "Nelder_Mead",
               optimizer = "nlminb",
               control = list(maxit = 1e6))

cat(paste("Length-weight parameter estimates for Eastern Bering pollock\n", 
          "Alpha:", exp(coef(wl_fit)[1]), "\n", "Beta:", 
          exp(coef(wl_fit)[2]), "\n", "sigma:", exp(coef(wl_fit)[3])))

pol.alpha <- as.numeric(exp(coef(wl_fit)[1]))
pol.beta <- as.numeric(exp(coef(wl_fit)[2]))

# 3.5- Plot model output versus observations
pol.dat$wl <- pred_wl(pol.dat$Length..mm., pol.alpha, pol.beta)

ggplot(pol.dat, aes(x = Length..mm.)) +
  geom_line(aes(y = wl), color = "red", size = 1.5) +
  geom_point(aes(y = Weight..gm.)) + 
  labs(title = "Observed (points) and predicted (line) weights at length for 
       EBS pollock", x = "Length (mm)", y = "Weight (gm)")

# 3.6- Fit LBV model to length-age data for EBS pollock
LVB_fit_pol <- mle2(NLL_LVB,
                    start = list(ln_Linf = log(max(pol.dat$Length..mm.)), 
                                 ln_k = log(0.2), t0 = 0, ln_sigma = log(0.5)),
                    data = list(obs.age = pol.dat$Age..years., 
                                obs.length = pol.dat$Length..mm.),
                    method = "Nelder-Mead",
                    optimizer = "nlminb",
                    control = list(maxit = 1e6))

cat(paste("LVB parameter estimates for Eastern Bering pollock\n", 
          "Linf estimate:", exp(coef(LVB_fit_pol)[1]), "\n", "k estimate:", 
          exp(coef(LVB_fit_pol)[2]), "\n", "t0 estimate:", coef(LVB_fit_pol)[3]), 
    "\n", "sigma estimate:", exp(coef(LVB_fit_pol)[4]))

Linf.pol <- as.numeric(exp(coef(LVB_fit_pol)[1]))
k.pol <- as.numeric(exp(coef(LVB_fit_pol)[2]))
t0.pol <- as.numeric(coef(LVB_fit_pol)[3])

# 3.7- Plot model output versus observations 
pol.dat$LVB <- pred_LVB(pol.dat$Age..years., Linf.pol, k.pol, t0.pol)

ggplot(pol.dat, aes(x = Age..years.)) +
  geom_line(aes(y = LVB), color = "red", size = 1.5) +
  geom_point(aes(y = Length..mm.)) + 
  labs(title = "Observed (points) and predicted (line) lengths at age for 
       EBS pollock", x = "Age", y = "Length (mm)")

# Problem 3- Continued...
  # Create vector of ages
ages <- c(1:15)

  # Predict weight at age for vector "ages"
waa <- pred_wl(pred_LVB(ages, Linf.pol, k.pol, t0.pol), pol.alpha, pol.beta)

  # Calculate expected number alive at each age
N <- vector(length = length(ages))

for(i in 1:(length(ages) - 1)) {
  N[1] <- 1000
  N[i + 1] <- N[i] * exp(-1 * 0.3)
}

  # Calculate expected biomass at age
B <- vector(length = length(ages))

for(i in 1:(length(ages))) {
  B[i] <- N[i] * waa[i]
}

  # Plot numbers, weight, and biomass at age
par(mar = c(4.1, 4.1, 1.5, 2.1),
    mfrow = c(3,1))
plot(N ~ ages, xlab = "", ylab = "Abundance", 
     main = "General natural mortality")
plot(waa ~ ages, xlab = "", ylab = "Weight (g)")
plot(B ~ ages, xlab = "Age", ylab = "Biomass (g)")

# Maximum biomass is reached at age 5 for this cohort

# Repeat above with different values for natural mortality 

  # Calculate expected number alive at each age
N.2 <- vector(length = length(ages))

for(i in 1:(length(ages) - 3)) {
  N.2[1] <- 1000
  N.2[2] <- N.2[1] * exp(-1 * 0.9)
  N.2[3] <- N.2[2] * exp(-1 * 0.45)
  N.2[i + 3] <- N.2[i + 2] * exp(-1 * 0.3)
}

# Calculate expected biomass at age
B.2 <- vector(length = length(ages))

for(i in 1:(length(ages))) {
  B.2[i] <- N.2[i] * waa[i]
}

# Plot numbers at age
par(mar = c(4.1, 4.1, 1.5, 2.1),
    mfrow = c(3,1))
plot(N.2 ~ ages, xlab = "", ylab = "Abundance", 
     main = "Age-specific natural mortality")
plot(waa ~ ages, xlab = "", ylab = "Weight (g)")
plot(B.2 ~ ages, xlab = "Age", ylab = "Biomass (g)")

# Maximum biomass is reached at age 6 for this cohort

