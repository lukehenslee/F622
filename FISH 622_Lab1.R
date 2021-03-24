#==================================================================================================
#Project Name: FISH 622 - QUANTITATIVE FISH POPULATION DYNAMICS: Lab #1
#Creator: Curry James Cunningham, College of Fisheries and Ocean Sciences, UAF
#Date: January 15, 2021
#
#Purpose: To practice fitting basic population dynamics models with R's function minimizers
#
# =========================================================================================

require(manipulate)
# install.packages("manipulate") # Install the manipulate package if you have not already done so.

library(bbmle)
# install.packages("bbmle")

# First, set your working directory!!!!!
setwd("/Users/curryc2/Documents/Students/2021/622 - Quant Fish Pop Dynamics/Content/Week 1/Lab")

# Part C: Simulating data from Myers and Worm (2003) ===============================

# First, we will see how I simulated a catch per unit effort time series similar to
#   those in Myers and Worm (2003)

# We can begin by defining a function in R that generates predictions, given values for 
#  the three parameters: N0, r, and delta.

# The arguments to this function will be these three parameters, plus a vector of time references (t)

pred_myers_worm <- function(N0=10, r=0.2, delta=0.5, time=c(0:50) ) {
  # First we generate our cpue predictions at each time point
  cpue <- N0*( (1-delta)*exp(-1*r*time) + delta )
  # As time is a vector, the resulting cpue object will be a vector of predicted values
  
  # Lets also have our function plot the predicted cpue time series
  plot(x=time, y=cpue, type="l", col="red", xlab="Time", ylab="cpue",
       ylim=c(0,max(cpue)))
  points(x=time, y=cpue, pch=21, bg="red")
  grid(col="black")
  
  # Finally, we can ask the function to return our predicted cpue values
  return(cpue)
}

# Lets test it out
time <- c(0:50)
time

# By calling the function we just created we plot the predicted cpue time series,
#   and return the predicted cpue values
cpue <- pred_myers_worm(N0=10, r=0.2, delta=0.1, time=time)
cpue

# Now, with some different parameter values
cpue <- pred_myers_worm(N0=10, r=0.05, delta=0.1, time=time)
cpue

# The manipulate function allows us to quickly create an interactive graphic
#   with sliders to control the parameters of our model


# CLICK THE LITTLE COG WHEEL IN THE UPPER-LEFT HAND CORNER OF YOUR PLOTTING WINDOW
#   TO MANIPULATE PARAMETER VALUES
library(manipulate)

manipulate(pred_myers_worm(N0, r, delta, time=time), 
           N0 = slider(min=1, max=25, initial=10, step=1),
           r = slider(min=0.001, max=1, initial=0.2, step=0.001), 
           delta = slider(min=0.001, max=0.999, initial=0.1, step=0.001)
         )

# Please explore how adjusting r and delta influence the shape of the relationship

# Now that we have a function to generate some predicted values, lets create
#   a time series with some random error.

# Fitting a model to data with no error isn't any fun!

# This can be done simply by generating a couple of cpue time series with different, and 
# adding random and normally-distributed errors
cpue_series1 <- pred_myers_worm(N0=10, r=0.2, delta=0.1, time=time)

set.seed(101) #I'll set the seed so we get the same results
cpue_series1_obs <- cpue_series1 + rnorm(n=length(time), mean=0, sd=0.3)

# How do these look?
plot(x=time, y=cpue_series1, type="l", xlab="Time", ylab="cpue", main="Simulated Data",
       ylim=c(0, max(cpue_series1, cpue_series1_obs)))
points(x=time, y=cpue_series1_obs, pch=21, bg=rgb(0,0,1, alpha=0.4))

# Lets generate another time series
cpue_series2 <- pred_myers_worm(N0=10, r=0.05, delta=0.2, time=time)
cpue_series2_obs <- cpue_series2 + rnorm(n=length(time), mean=0, sd=0.3)
plot(x=time, y=cpue_series2, type="l", xlab="Time", ylab="cpue", main="Simulated Data",
     ylim=c(0, max(cpue_series2, cpue_series2_obs)))
points(x=time, y=cpue_series2_obs, pch=21, bg=rgb(0,0,1, alpha=0.4))

# An we will package these and save them as a .csv
sim.cpue <- data.frame(time, cpue_series1_obs, cpue_series2_obs)
write.csv(sim.cpue, "sim.cpue.csv")

# Now, simulate your own data, but increase the sd= in the rnorm() call when
#   you add random errors. Then copy the time series into your spreadsheet and 
#   try fitting the model again.
# What happens to your estimated parameter values, relative to their "true" values?

# Part D: Fitting the Modified Exponential in R =================================

# Now that we have seen how we can fit the modified exponential model to simulated
#   cpue time series in Excel, lets do it in R!

# Much like Solver in Excel, the mle2() function in R offers an efficient
#   non-linear function minimizer. 

library(bbmle)

?mle2

# For this we will need both our simulated cpue data, and a new function that both
#   1) Generates predicted values
#   2) Calculates the sum of squared deviations between predicted and observed values

# The arguments to this function will be values for the parameters of the model,
#   our time references, and our observed cpue time series (obs.cpue)

ssq_myers_worm <- function(N0=10, r=0.2, delta=0.1, time=time, obs.cpue=NULL) {
  # First, generate predictions based on parameter values
  pred.cpue <- cpue <- N0*( (1-delta)*exp(-1*r*time) + delta )
  
  # Calculate squared deviations between observed and predicted values, at each time point
  sq <- (obs.cpue-pred.cpue)^2
  
  # Calculate and return the sum of squared deviations
  ssq <- sum(sq)
  return(ssq)
}

# Ok, lets test it out... 

# We can simulate data and add error as before, with known values for the parameters

time <- c(0:50)

# REMEMBER THESE "TRUE" VALUES FOR N0, r, and delta.
#   We will compare these with those estimated by the function minimizer a little 
#   further down. 
obs.cpue <- pred_myers_worm(N0=15, r=0.2, delta=0.1, time=time) + rnorm(n=length(time), mean=0, sd=0.5)

plot(x=time, y=obs.cpue, type="p", xlab="Time", ylab="cpue", main="Simulated Data",
     ylim=c(0, max(obs.cpue)), pch=21, bg=rgb(0,0,1, alpha=0.3))

# Test our SSQ function
ssq_myers_worm(N0=15, r=0.2, delta=0.1, time=time, obs.cpue=obs.cpue)

# What happens if we move our input r value away from the "true" value
#   that was used to simulate obs.cpue?
ssq_myers_worm(N0=15, r=0.5, delta=0.1, time=time, obs.cpue=obs.cpue)
ssq_myers_worm(N0=15, r=0.1, delta=0.1, time=time, obs.cpue=obs.cpue)

# Perfect, it gets higher, that is what we want. 

# Now, to actually do the minimization function we need to provide several things
#   to the mle2() function

# The first argument is our function that returns a value to be minimized 
#   (the objective value)
# The next critical argument is the list of starting values for the parameters.
#   This also tells mle2() which of the arguments to our function it should it should search across
#  Finally, we specify a list indicating which inputs to our function to be minimzed
#   are data

fit <- mle2(ssq_myers_worm, 
            start=list(N0=12, r=0.1, delta=0.5),
            data=list(time=time, obs.cpue=obs.cpue),
            method="Nelder-Mead")

fit

# Coefficients as a named vector
coef(fit)

# Nice! It looks like it was able to estimate or parameter values with reasonable precision.

# Part E: Fitting the Logistic Model to US Population Data in R =============================

# First we will load the US Population data 

us.dat <- read.csv("US Population.csv", header=TRUE)

str(us.dat)

# Lets visualize the US population data
plot(humans~year, data=us.dat, type="h", lwd=2, col="blue",
       main="US Population", xlab="Year", ylab="Human (millions)")
grid()
points(humans~year, data=us.dat, pch=21, bg="blue")

# Next, lets create a function to generate predicted values from the logistic model

logistic_fxn <- function(N0, r, K, years) {
  pred <- (N0*exp(r*(years-1790))) / (1-(N0/K)+(N0/K)*exp(r*(years-1790)))
  return(pred)
}

# Does it work?
logistic_fxn(N0=10, r=0.1, K=500, years=us.dat$year)

# Now we will wrap that in a function to calculate the SSQ
logistic_ssq <- function(N0, r, K, years, obs.humans) {
  # Predicted values (calling the function above)
  pred <- logistic_fxn(N0=N0, r=r, K=K, years=years)
  
  # Squared deviations
  sq <- (obs.humans - pred)^2
  
  # Sum of squared deviations
  ssq <- sum(sq)
  return(ssq)
}

# We can make sure this works...
logistic_ssq(N0=10, r=0.1, K=500, years=us.dat$year, obs.humans=us.dat$humans)
logistic_ssq(N0=10, r=0.2, K=500, years=us.dat$year, obs.humans=us.dat$humans) # Worse
logistic_ssq(N0=10, r=0.3, K=500, years=us.dat$year, obs.humans=us.dat$humans) # Worse
logistic_ssq(N0=10, r=0.4, K=500, years=us.dat$year, obs.humans=us.dat$humans) # Worse

logistic_ssq(N0=10, r=0.05, K=500, years=us.dat$year, obs.humans=us.dat$humans) # Better
logistic_ssq(N0=10, r=0.01, K=500, years=us.dat$year, obs.humans=us.dat$humans) # Worse

# From this we can see that r values close to 0.05 provide better fits, compared to r=0.01 or 0.1

# While it is not good practice to fit models by eye, this can give us a good idea
#   of reasonable starting values
# Observed
plot(humans~year, data=us.dat, type="h", lwd=2, col="blue",
     main="US Population", xlab="Year", ylab="Human (millions)")
grid()
points(humans~year, data=us.dat, pch=21, bg="blue")
# Predicted
pred <- logistic_fxn(N0=10, r=0.1, K=500, years=us.dat$year)
lines(x=us.dat$year, y=pred, col="red")

# NOT GREAT!

# Observed
plot(humans~year, data=us.dat, type="h", lwd=2, col="blue",
     main="US Population", xlab="Year", ylab="Human (millions)")
grid()
points(humans~year, data=us.dat, pch=21, bg="blue")
# Predicted
pred <- logistic_fxn(N0=10, r=0.05, K=500, years=us.dat$year)
lines(x=us.dat$year, y=pred, col="red")

# BETTER

# Observed
plot(humans~year, data=us.dat, type="h", lwd=2, col="blue",
     main="US Population", xlab="Year", ylab="Human (millions)")
grid()
points(humans~year, data=us.dat, pch=21, bg="blue")
# Predicted
pred <- logistic_fxn(N0=10, r=0.01, K=500, years=us.dat$year)
lines(x=us.dat$year, y=pred, col="red")

# TOO LOW!

# Observed
plot(humans~year, data=us.dat, type="h", lwd=2, col="blue",
     main="US Population", xlab="Year", ylab="Human (millions)")
grid()
points(humans~year, data=us.dat, pch=21, bg="blue")
# Predicted
pred <- logistic_fxn(N0=10, r=0.02, K=500, years=us.dat$year)
lines(x=us.dat$year, y=pred, col="red")

# NOT BAD

# OK, NOW THAT WE HAVE SOME DECENT STARTING VALUES, LETS FIT THE MODEL

fit <- mle2(logistic_ssq,
            start=list(N0=10, r=0.02, K=500),
            data=list(years=us.dat$year, obs.humans=us.dat$humans),
            method="Nelder-Mead")

fit
coef(fit)



# How to these compare to what Solver estimated in Part A?

# Lets plot the model fit:
# Observed
plot(humans~year, data=us.dat, type="h", lwd=2, col="blue",
     main="US Population", xlab="Year", ylab="Human (millions)")
grid()
points(humans~year, data=us.dat, pch=21, bg="blue")
# Predicted
pred <- logistic_fxn(N0=coef(fit)["N0"], r=coef(fit)["r"], K=coef(fit)["K"], years=us.dat$year)
lines(x=us.dat$year, y=pred, col="red")


# Turns out even with rather poor starting values, we still arrive at the same solution
# mle2() is a very robust function minimizer for a simple problem like this!
fit <- mle2(logistic_ssq,
            start=list(N0=20, r=0.1, K=200),
            data=list(years=us.dat$year, obs.humans=us.dat$humans),
            method="Nelder-Mead")

fit
coef(fit)



