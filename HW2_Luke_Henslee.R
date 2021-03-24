#==================================================================================================
# FISH 622 homework #2: Surplus production
# 3/1/2021
# Luke Henslee, College of Fisheries and Ocean Sciences, UAF
#
#==================================================================================================
#NOTES:

#==================================================================================================

# Load packages and set working directory ####
library(tidyverse)
library(ggplot2)
library(bbmle)

setwd("C:/Users/lhhenslee/Desktop/Luke/School/F622")
hw2.data <- file.path(getwd(), "data")
hw2.data


# Problem 1 ####

# 1.1- Load dataset and visualize
lat_prod_dat <- read.csv(file=file.path(hw2.data, "Latent Productivity Data.csv"))

attach(lat_prod_dat)
plot(Biomass, Pdot)

# 1.2- Set maxP
maxP <- max(Pdot)
  # A reasonable value for Binf might be around 110

# 1.3- Create prediction functions
  # 1.3a- Graham-Schaefer
sim_Pdot_GS <- function(B, m, Binf) {
  P <- (4*m/Binf) * (1-(B/Binf)) * B
  return(P)
}

  # 1.3b- Pella-Tomlinson
sim_Pdot_PT <- function (B, m, Binf, n) {
  gamma <- (n^(n/(n-1))) / (n-1)
  P <- gamma * m * (B/Binf) - gamma * m * ((B/Binf)^n)
  return(P)
}

  # 1.3c- Fox
sim_Pdot_Fox <- function(B, m, Binf) {
  P <- -exp(1) * m * (B/Binf) * log(B/Binf)
  return(P)
}

# 1.4- Create negative log-likelihood functions for each model
  # 1.4a- Graham-Schaefer
NLL_GS <- function(ln_m, ln_Binf, B, obs_Pdot, ln_sigma) {
  # Exponentiate model parameters
  m <- exp(ln_m)
  Binf <- exp(ln_Binf)
  sigma <- exp(ln_sigma)
  # Calculate predictions under GS model
  pred_Pdot <- sim_Pdot_GS(B, m, Binf)
  # Calculate log likelihood using dnorm()
  NLL <- -1 * dnorm(obs_Pdot, pred_Pdot, sigma, log = TRUE)
  return(sum(NLL))
}

  # 1.4b- Pella-Tomlinson
NLL_PT <- function(ln_m, ln_Binf, ln_n, B, obs_Pdot, ln_sigma) {
  # Exponentiate model parameters
  m <- exp(ln_m)
  Binf <- exp(ln_Binf)
  sigma <- exp(ln_sigma)
  n <- exp(ln_n)
  # Calculate predictions under PT model
  pred_Pdot <- sim_Pdot_PT(B, m, Binf, n)
  # Calculate log likelihood using dnorm()
  NLL <- -1 * dnorm(obs_Pdot, pred_Pdot, sigma, log = TRUE)
  return(sum(NLL))
}

  #1.4c- Fox
NLL_Fox <- function(ln_m, ln_Binf, B, obs_Pdot, ln_sigma) {
  # Exponentiate model parameters
  m <- exp(ln_m)
  Binf <- exp(ln_Binf)
  sigma <- exp(ln_sigma)
  # Calculate predictions under Fox model
  pred_Pdot <- sim_Pdot_Fox(B, m, Binf)
  # Calculate log likelihood using dnorm()
  NLL <- -1 * dnorm(obs_Pdot, pred_Pdot, sigma, log = TRUE)
  return(sum(NLL))
}

# 1.6- Use mle2() to fit each model
  # 1.6a- Graham-Schaefer
GS_fit <- mle2(NLL_GS,
               start = list(ln_m = log(maxP), ln_Binf = log(120), ln_sigma = log(2)),
               data = list(B = lat_prod_dat$Biomass, obs_Pdot = lat_prod_dat$Pdot),
               method = "Nelder-Mead")
  # Save coefficients
GS_coef <- exp(coef(GS_fit))
GS_m <- as.numeric(exp(coef(GS_fit)[1]))
GS_Binf <- as.numeric(exp(coef(GS_fit)[2]))


  # 1.6b- Pella_Tomlinson
PT_fit <- mle2(NLL_PT,
               start = list(ln_m = log(maxP), ln_Binf = log(120), ln_sigma = log(2), ln_n = log(4)),
               data = list(B = lat_prod_dat$Biomass, obs_Pdot = lat_prod_dat$Pdot),
               method = "Nelder-Mead")
  # Save coefficients
PT_coef <- exp(coef(PT_fit))
PT_m <- as.numeric(exp(coef(PT_fit)[1]))
PT_Binf <- as.numeric(exp(coef(PT_fit)[2]))
PT_n <- as.numeric(exp(coef(PT_fit)[3]))

  # 1.6c- Fox
Fox_fit <- mle2(NLL_Fox,
               start = list(ln_m = log(maxP), ln_Binf = log(120), ln_sigma = log(2)),
               data = list(B = lat_prod_dat$Biomass, obs_Pdot = lat_prod_dat$Pdot),
               method = "Nelder-Mead")
  # Save coefficients
Fox_coef <- exp(coef(Fox_fit))
Fox_m <- as.numeric(exp(coef(Fox_fit)[1]))
Fox_Binf <- as.numeric(exp(coef(Fox_fit)[2]))

# 1.7 Plot data as points and model predictions as lines
  # Create columns of model predicted values for Pdot
lat_prod_dat$GS <- sim_Pdot_GS(Biomass, GS_m, GS_Binf)

lat_prod_dat$PT <- sim_Pdot_PT(Biomass, PT_m, PT_Binf, PT_n)

lat_prod_dat$Fox <- sim_Pdot_Fox(Biomass, Fox_m, Fox_Binf)
  
  # Plot
ggplot(lat_prod_dat, aes(x = Biomass)) +
  geom_point(aes(y = Pdot)) +
  geom_line(aes(y = GS, color = "GS"), size = 1.5) +
  geom_line(aes(y = PT, color = "PT"), size = 1.5) +
  geom_line(aes(y = Fox, color = "Fox"), size = 1.5, linetype = "dashed") +
  labs(x = "Biomass", y = "Pdot", color = "Model fit") +
  scale_color_manual(values = c("GS" = "blue", "PT" = "red", "Fox" = "green"))

  # From these plots, it seems that the Fox and Pella-Tomlinson models are about
  # the same, and both fit the data better than the Graham-Shaefer model. 

# 1.8- Compare and contrast parameter estimates
GS_coef
PT_coef
Fox_coef
  # The first thing I notice is that the m coefficient is almost the same in the 
  # PT and Fox models, but the estimate for m in the GS model is a closer fit to
  # our data. I visually estimated Binf to be about 110 from our data values,
  # and all models were about 6 units off from that. Sigma estimates are about the 
  # same for the PT and Fox models, and much higher for the GS model. 

# 1.9- Use AIC() for each model
AIC(GS_fit, PT_fit, Fox_fit)
  # The Fox model has the lowest AIC score, but it is very close to the PT model
  # score. The GS model is way off, and scores about 150 points higher than 
  # either of the other two models. The Fox model provides the most parsimonious 
  # fit to our data. Fox has fewer parameters...

# Problem 2 ####
# Load data and visualize
hake_dat <- read.table(file=file.path(hw2.data, "hake.dat"), header = TRUE)

attach(hake_dat)
par(mar = c(5, 4, 4, 4) + 0.1)
plot(Year, Catch, type = "l", col = "green")
par(new = TRUE)
plot(Year, Index, type = "l", col = "red", axes = F, ylab = "")
mtext(" Biomass index", side = 4, line = 3)
axis(4)
legend("topright", legend=c("Catch", "Biomass index"),
       col=c("green", "red"), lty = 1, cex = 0.8)
  # It seems that in the beginning of the time series, increased catch leads to 
  # biomass index decreases until the mid 1970's, when catch begins falling along
  # with index. Around 1980, the catch and biomass index start to look more 
  # directly related instead of having a negative correlation. 

# 2.1- Create predictive function
biom_logistic <- function(Ct, n.years, r, K) {
  B <- vector(length = n.years)
  B[1] <- K
  for(i in (1:(n.years-1))) {
    B[i+1] <- B[i] + (r*B[i] * (1 - (B[i]/K))) - Ct[i]
    B[i+1] = max(B[i+1], 1e-3)
  }
  return(B)
}

# 2.2- Create SSQ function
ssq_logistic <- function(ln_q, ln_r, ln_K, n.years, Ct, obsI) {
  q <- exp(ln_q)
  r <- exp(ln_r)
  K <- exp(ln_K)
  predI <- q*biom_logistic(Ct, n.years, r, K)
  SSQ <- sum((log(obsI) - log(predI))^2)
  return(SSQ)
}

# 2.3- Use mle2 to estimate parameter values
biom_fit <- mle2(ssq_logistic,
                 start = list(ln_q = log(1e-5), ln_r = log(0.2), ln_K = log(2000)),
                 data = list(n.years = length(hake_dat$Year), Ct = hake_dat$Catch, 
                             obsI = hake_dat$Index),
                 method = "Nelder-Mead")

  # Save parameter values
Log_I_coef <- exp(coef(biom_fit))
Log_q <- Log_I_coef[1]
Log_r <- Log_I_coef[2]
Log_K <- Log_I_coef[3]

# 2.4 Plot data as points and model predictions as lines  

  # Create column of model predicted values for Ihat
hake_dat$Log_I_pred <- Log_q * biom_logistic(hake_dat$Catch, length(hake_dat$Year), Log_r, Log_K)

  # Plot
ggplot(hake_dat, aes(x = Year)) +
  geom_point(aes(y = Index)) +
  geom_line(aes(y = Log_I_pred, color = "Log_I_pred"), size = 1) +
  labs(x = "Year", y = "Biomass index", color = "Model fit") +
  scale_color_manual(values = c("Log_I_pred" = "blue"))

# 2.5 Plot surplus production as a function of biomass
  # Create function to simulate surplus production and plot
Log_surplus <- function(r, K) {
  # A reasonable range of biomass would be from 0 to our model estimate of carrying
  # capacity, then a little further for good measure... 
  Bt <- as.vector(c(0:(K + 200)))
  SP <- vector(length = length(Bt))
  for(i in (1:length(Bt))) {
    SP[i] <- r*Bt[i] * (1 - (Bt[i]/K))
  }
  plot(Bt, SP, type = "l")
  abline(h = max(SP), lty = 2, col = "red")
  abline(v = K/2, lty = 2, col = "blue")
  legend("topright", legend=c("MSY", "Bmsy"),
         col=c("red", "blue"), lty = 2, cex = 0.8)
}

  # Generate simulated data from model parameter estimates of r and K
Log_surplus(Log_r, Log_K)

  # Visually, it appears that MSY is about 250 and the biomass producing MSY is
  # approximately 1400

# 2.6- Use analytical solution to find MSY and Bmsy
Log_MSY_Bmsy <- function(r, K) {
  MSY <- (r*K)/4
  Bmsy <- K/2
  return(as.numeric(c(MSY, Bmsy)))
}

Log_MSY_Bmsy(Log_r, Log_K)

  # Our analytical value for MSY is 261.32 and Bmsy is 1416.67

# Problem 2b- Repeat with Fox model ####
# 2.1b- Create predictive function
biom_fox <- function(Ct, n.years, r, K) {
  B <- vector(length = n.years)
  B[1] <- K
  for(i in (1:(n.years-1))) {
    B[i+1] <- B[i] + (r*B[i] * (1 - (log(B[i])/log(K)))) - Ct[i]
    B[i+1] = max(B[i+1], 1e-3)
  }
  return(B)
}

# 2.2b- Create SSQ function
ssq_fox <- function(ln_q, ln_r, ln_K, n.years, Ct, obsI) {
  q <- exp(ln_q)
  r <- exp(ln_r)
  K <- exp(ln_K)
  predI <- q*biom_fox(Ct, n.years, r, K)
  SSQ <- sum((log(obsI) - log(predI))^2)
  return(SSQ)
}

# 2.3b- Use mle2 to estimate parameter values
biom_fit_fox <- mle2(ssq_fox,
                 start = list(ln_q = log(5e-4), ln_r = log(1.7), ln_K = log(5000)),
                 data = list(n.years = length(hake_dat$Year), Ct = hake_dat$Catch, 
                             obsI = hake_dat$Index),
                 method = "Nelder-Mead")

# Save parameter values
Fox_I_coef <- exp(coef(biom_fit_fox))
Fox_q <- Fox_I_coef[1]
Fox_r <- Fox_I_coef[2]
Fox_K <- Fox_I_coef[3]

# 2.4b Plot data as points and model predictions as lines  

# Create column of model predicted values for Ihat
hake_dat$Fox_I_pred <- Fox_q * biom_fox(hake_dat$Catch, length(hake_dat$Year), Fox_r, Fox_K)

# Plot
ggplot(hake_dat, aes(x = Year)) +
  geom_point(aes(y = Index)) +
  geom_line(aes(y = Fox_I_pred, color = "Fox_I_pred"), size = 1) +
  labs(x = "Year", y = "Biomass index", color = "Model fit") +
  scale_color_manual(values = c("Fox_I_pred" = "red"))

# 2.5b Plot surplus production as a function of biomass
# Create function to simulate surplus production and plot
Fox_surplus <- function(r, K) {
  # A reasonable range of biomass would be from 1 to our model estimate of carrying
  # capacity, then a little further for good measure... 
  Bt <- as.vector(c(1:(K + 200)))
  SP <- vector(length = length(Bt))
  for(i in (1:length(Bt))) {
    SP[i] <- r*Bt[i] * (1 - (log(Bt[i])/log(K)))
  }
  plot(Bt, SP, type = "l")
  abline(h = max(SP), lty = 2, col = "red")
  abline(v = 0.37*K, lty = 2, col = "blue")
  legend("topright", legend=c("MSY", "Bmsy"),
         col=c("red", "blue"), lty = 2, cex = 0.8)
}

# Generate simulated data from model parameter estimates of r and K
Fox_surplus(Fox_r, Fox_K)

# Visually, it appears that MSY is about 250 and the biomass producing MSY is
# approximately 1200

# 2.6b- Use analytical solution to find MSY and Bmsy
Fox_MSY_Bmsy <- function(r, K) {
  # Peak catch occurs when B = 0.37K
  Bmsy <- 0.37 * Fox_K
  # Plug Bmsy in for Bt in the Fox surplus production model 
  MSY <- r * Bmsy * (1 - log(Bmsy) / log(K))
  return(as.numeric(c(MSY, Bmsy)))
}

Fox_MSY_Bmsy(Fox_r, Fox_K)

# Time allocation ####

  # I would estimate this assignment took me about six hours, including a one-
# hour study session with some classmates. 