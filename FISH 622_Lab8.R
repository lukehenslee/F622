#==================================================================================================
#Project Name: FISH 622 - QUANTITATIVE FISH POPULATION DYNAMICS: Lab #8 YPR and Matrix Models
#Creator: Curry James Cunningham, College of Fisheries and Ocean Sciences, UAF
#Date: March 26, 2021
#
#Purpose: 
#
# =========================================================================================
library(tidyverse)
library(dplyr)
library(manipulate)
library(bbmle)
library(visreg) #Install any of these packages if you don't have them already
library(ggthemes)
library(manipulate)
library(ggridges)
library(reshape2)
library(matrixcalc)

# Set your working directory!!!!!
setwd("/Users/curryc2/Documents/Students/2021/622 - Quant Fish Pop Dynamics/Content/Week 8/Lab")

# Part A: Tracking a Cohort ==================================================

# In this exercise, we will repeat the process we followed in Lab #7 to track the decline
#  in the abundance of a cohort across subsequent ages.

# But, we will show how to do so in R, using loops to iterate across ages
#   and if/else (logic) statements 

# First we need define the ages of key life history events for a commercially important
#   species (Muppetus muppetus)

# The assumed age at recruitment
tr <- 2 # years

# The age at which fish become available to the fishery
tc <- 10 # years

# One assumption inherent in the Beverton-Holt approximation is "knife-edged"
#   selectivity, where there is no exploitation of the cohort
#   at ages below tc, and full exploitation at all ages above tc.

# We also need to define our initial number of recruits for the cohort.
#   This is somewhat arbitrary

# Here we will define it as 1,000 individuals for the cohort alive at the time of recruitment
Nr <- 1000

# Finally, we need the natural mortality rate for M. muppetus
M <- 0.1

# And an assumed fishing mortality
Fmort <- 0.1

# From this we can calculate total instantaneous mortalty: Z=F+M
Z <- Fmort + M

# Given values for Nr, M, tc and tr, we can calculate the abundance of the cohort 
#   at when they first become available to the fishery
Nc <- Nr*exp(-M*(tc-tr))
Nc

# We can see that the survival fraction from recruitment to tc:
Nc/Nr

# To show this, we can plot out the decline in our cohort as we did before (Lab #7),
#   across a range of ages

trial.ages <- 1:50
trial.ages

# Get length of vector of trial ages
n.trial.ages <- length(trial.ages)

# Lets create a vector to hold our predicted abundances across ages for our cohort
Npred <- vector(length=n.trial.ages)

# Lets loop over ages and calculate our abundance at age
#   We will use a for() loop to do so, but will need to embed an if()
#   statement to account for the impact of M vs. Z=F+M, before/after tc

# In this for loop, our iterator "age" will take the value of each
#   element in our sequence "trial.ages", and repeat a number of times equal to
#   the length of our sequence.

# Just for clarity, lets prove this to ourselves

i <- 1
for(t in trial.ages) {
  print(paste("age:", t))
  print(paste("Number of times our loop has executed:", i))
  i <- i+1
}

# OK, now that we hopefully have a better sense of how the loop will execute
#   we can use this same loop to populate our expected numbers-at-age for the cohort.
# REMEMBER: This is a lot like dragging our cell reference down across ages in Excel!

for(t in trial.ages) {
  if(t<tc) {
    Npred[t] <- Nr*exp(-M*(t-tr))
  }else {
    Npred[t] <- Nr*exp(-M*(tc-tr))*exp(-Z*(t-tc))
  }
}

# Lets plot out what we calculated
plot(x=trial.ages, y=Npred, type="l", col="blue",
       xlab="Age", ylab="Cohort Abundance: N(t)", lwd=2)

# Age at recruitment
segments(x0=tr, y0=-100, x1=tr, y1=Nr, col="red", lty=3)
segments(x0=-100, y0=Nr, x1=tr, y1=Nr, col="red", lty=3)
points(x=tr, y=Nr, pch=21, bg="red")

# Age at exploitation (fishery)
segments(x0=tc, y0=-100, x1=tc, y1=Nc, col="orange", lty=3)
segments(x0=-100, y0=Nc, x1=tc, y1=Nc, col="orange", lty=3)
points(x=tc, y=Nc, pch=21, bg="orange")

# Next, lets approximate the first derivative of this N(t) curve with respect to age t
#   i.e. the slope of the line at each time point

# To approximate dN/dt we can substract subsequent abundances from the prior abundance
dN.dt <- vector(length=n.trial.ages)
dN.dt

for(t in 1:(n.trial.ages-1)) {
  dN.dt[t] <- Npred[t+1]-Npred[t]
}

dN.dt
# This vector should be equal to dN/dt should be equal to -M*N for ages less than tc
#   and equal to -Z*M for t >= tc.

# Or equivalently the per-capita rate of change in abundance dN/Ndt should be equal to -M before tc, and equal to -Z 

dN.dt/Npred

# Lets plot it...
plot(x=trial.ages, y=dN.dt/Npred, type='b')

# Remember this is an approximation to a continuous process, so it isn't exact.
#   WHY?

Z
M

tc

# Part B: Beverton-Holt LVB Isometric Yield Calculations =========================================

# Beverton and Holt (1957) showed that yield can be calculated analytically
#   in a (somewhat) straight froward fashion if the weight-age relationship takes the isometric
#   von Bertalanffy (LVB) form. 

# In general terms this means that we can track a cohort across its fishable lifespan
#   and determine the total yield for the cohort, based on co-occurring
#   mortality and growth processes (i.e. abundance of the cohort declining across ages
#   and the average weight of individuals within the cohort increasing across ages).


# Lets begin by defining some parameters of the weight-age relationship for a 
#  commercially important species (Muppetus muppetus)
Winf <- 2 # in Kg
k <- 0.2
t0 <- 0

# Next, lets plot out our assumed weight-age relationship, across a range of ages
ages <- 0:50

# Expected weight at age
waa <- Winf*(1-exp(-k*(ages-t0)))^3

# Plot weight-age relationship
plot(x=ages, y=waa, type="l", col="red",
     xlab="Age (years)", ylab="Average Weight (kg)")
grid()
points(x=ages, y=waa, pch=21, bg="red")

# In order to calculate YPR we also need to define several key parameters

# Using a cubic expansion for the weight-age relationship we can define total (cumulative)
#   yield (in units of biomass i.e. kg) from a cohort
#   up through a given age (t)

# Lets create a function to calculate total yield through a given age (t)
#   for the cohort.

# We will break this calculation into steps within the function

BH_yield <- function(t, Fmort, M,
                     Winf, k, t0,
                     tr, tc) {

  # Lets calculate total mortality
  Z <- Fmort + M
  
  # The parameter Un for the cubic expansion takes the following values across 
  #   the summation 0-3
  Un <- c(1,-3,3,-1)
  
  # Finally, lets calculate Nc, the cohort abundance at the age it becomes 
  #   available to the fishery
  Nc <- Nr*exp(-M*(tc-tr))
  
  # For the total yield at age t calculations, we can break it into two parts
  #   first the component before the summartion
  Yt.part1 <- Fmort*Nc*Winf
  
  # Next, we need to do the summation ...
  Yt.part2 <- 0
  
  # We will use a "counter variable" i to reference locations in Un
  i <- 1
  for(n in c(0:3)) {
    Yt.part2 <- Yt.part2 + ( Un[i] / (Z+n*k) ) * 
                             exp(-n*k*(tc-t0)) * 
                             (1-exp(-(Z+n*k) * (t-tc) ))
    i <- i+1 #Update the counter variable
  } #next n

  #The total yield up to age t is then
  Yt <- Yt.part1 * Yt.part2
  
  # Return the total Yield
  return(Yt)
}
  
# Lets test our function
tc

BH_yield(t=3, Fmort=Fmort, M=M,
         Winf=Winf, k=k, t0=t0,
         tr=tr, tc=tc)

# Oh no we specified an age (t) less than tc.. we can't calculate yield as there is none

# At the start of age tc, we have no catch
BH_yield(t=10, Fmort=Fmort, M=M,
         Winf=Winf, k=k, t0=t0,
         tr=tr, tc=tc)

# But we can calculate total yield Y(t) for subsequent ages

BH_yield(t=11, Fmort=Fmort, M=M,
         Winf=Winf, k=k, t0=t0,
         tr=tr, tc=tc)
BH_yield(t=15, Fmort=Fmort, M=M,
         Winf=Winf, k=k, t0=t0,
         tr=tr, tc=tc)
BH_yield(t=50, Fmort=Fmort, M=M,
         Winf=Winf, k=k, t0=t0,
         tr=tr, tc=tc)

# Lets simulate total yield across a range of ages

trial.yield.ages <- c(tc:50)
trial.yield.ages

yield <- BH_yield(t=trial.yield.ages, Fmort=Fmort, M=M,
                  Winf=Winf, k=k, t0=t0,
                  tr=tr, tc=tc)

# Lets plot this out
plot(x=trial.yield.ages, y=yield, type="l", col="red", lwd=2,
     ylab="Total Yield to Age t: Y(t)",
     xlab="Age (t)")
points(x=trial.yield.ages, y=yield, pch=21, bg="red")
grid()

#Lets plot cohort abundance and yield at age
plot(x=trial.ages, y=Npred, type="l", col="blue", lwd=3,
       xlab="Age (t)", ylab="Abundance or Yield")
points(x=trial.ages, y=Npred, pch=21, bg="blue")

# Total Yield
lines(x=trial.yield.ages, y=yield, col="red", lwd=2)
points(x=trial.yield.ages, y=yield, pch=21, bg="red")
legend("topright", legend=c("Abundance: N(t)", "Total Yield: Y(t)"), text.col=c("blue","red"))

# Next, lets create a "wrapper function" that generates predictions and 
#   creates the plot, so we can feed it into manipulate()

plot_BH_yield <- function(trial.ages=c(1:50), 
                          Fmort=0.1, M=0.1,
                          Winf=2, k=0.2, t0=0,
                          tr=2, tc=10) {
  par(mfrow=c(2,1), oma=c(0,0,0,0), mar=c(4,4,1,1))
  # Plot Weight at Age
  waa <- Winf*(1-exp(-k*(trial.ages-t0)))^3
  plot(x=trial.ages, y=waa, type="l", col="red",
       xlab="Age (t)", ylab="Average Weight (kg)")
  grid()
  points(x=trial.ages, y=waa, pch=21, bg="red")
  
  # Plot Abundance at Age and Total Yield at Age
  
  # Calculate Abundance at age: N(t)
  Npred <- vector(length=length(trial.ages))
  
  # Calculate total instantaneous mortality
  Z <- Fmort + M
  
  for(t in trial.ages) {
    if(t<tc) {
      Npred[t] <- Nr*exp(-M*(t-tr))
    }else {
      Npred[t] <- Nr*exp(-M*(tc-tr))*exp(-Z*(t-tc))
    }
  }
  
  # Calculate total (cumulative) yield at age: Y(t)
  trial.yield.ages <- tc:max(trial.ages)
  yield <- BH_yield(t=trial.yield.ages, Fmort=Fmort, M=M,
                    Winf=Winf, k=k, t0=t0,
                    tr=tr, tc=tc)
  
  # Plotting N(t) and Y(t)
  plot(x=trial.ages, y=Npred, type="l", col="blue", lwd=3,
       xlab="Age (t)", ylab="Abundance or Yield")
  points(x=trial.ages, y=Npred, pch=21, bg="blue")
  
  # Total Yield
  lines(x=trial.yield.ages, y=yield, col="red", lwd=2)
  points(x=trial.yield.ages, y=yield, pch=21, bg="red")
  legend("topright", legend=c("Abundance: N(t)", "Total Yield: Y(t)"), text.col=c("blue","red"))
}

# First, lets make sure our function works
dev.off()

plot_BH_yield(trial.ages=c(1:50), 
                          Fmort=0.1, M=0.1,
                          Winf=2, k=0.2, t0=0,
                          tr=2, tc=10)

plot_BH_yield(trial.ages=c(1:50), 
              Fmort=0.1, M=0.1,
              Winf=2, k=0.2, t0=0,
              tr=2, tc=5)

# OK, lets use manipulate to explore how Abundance and Yield change
#   As a function of M, F, and our growth parameters

manipulate(plot_BH_yield(trial.ages=c(1:50), 
                         Fmort, M,
                         Winf, k, t0,
                         tr=2, tc),
           Fmort=slider(min=0, max=0.5, initial=0.1, step=0.001),
           M=slider(min=0, max=0.5, initial=0.1, step=0.001),
           
           tc=slider(min=tr, max=50, initial=5, step=1),
           
           Winf=slider(min=0, max=5, initial=2, step=0.001),
           k=slider(min=0, max=1, initial=0.2, step=0.001),
           t0=slider(min=-5, max=5, initial=0.1, step=0.001)
           )


# Part C: Beverton-Holt LVB Isometric YPR Calculations =========================================

# In Part B we explored how we can use the Beverton-Holt (1957) approximation to 
#   simulate total yield from a cohort through a given age Y(t).

# Next we will calculate yield-per-recruit, a reference metric that is useful
#   when considering the influence of the length at which fish become 
#   available to harvest (Lc) and the effective instantaneous fishing mortality
#   rate (F)

# This formulation in terms of a minimum length limit (Lc) instead of a minimum
#   age limit (tc) is convenient, because as a harvester you are unlikely
#   to know the age of a given fish, but you can easily measure it, and either retain
#   or return it to the water. NOTE: We will ignore post-release (discard) 
#   mortality ... for now.

# Beverton-Holt (1957) showed that lifetime yield from a cohort (y) can be 
#   approximated assuming an isometric weight-age relationship
# As a function of the "exploitation rate" (E=F/Z) and several derived parameters:

# 1) The ratio of the minimum size limit to the asymptotic maximum length (c = Lc/Linf)
# 2) The ratio of natural mortality to the (k) growth parameter (m = M/k)

# Lets create a function to calculate yield per recruit.

BH_ypr <- function(Fmort, M, tr,
                   c=0.5, # Ratio of Lc/Linf
                   Winf, k, t0) {
  # First, lets calculate total instantaneous
  Z <- Fmort + M
  
  # Second, lets calculate the exploitation rate
  E <- Fmort/Z
  
  # Third, we will define the derived variable (m) as the ratio of natural mortality (M)
  #   to the growth coefficient (k)
  m <- M/k
  
  # Fourth, define parameter Un for the cubic expansion takes the following values across 
  #   the summation 0-3
  Un <- c(1,-3,3,-1)
  
  # Fifth, we will calculate "lifetime yield" from the cohort
  #   As before we need to do a summation across values for n for the cubic expansion
  
  # We will define a temporary variable whose value we will update (by adding)
  #   in each step of our loop to do the summarion.
  y.temp <- 0
  
  # We will use a "counter variable" i to reference locations in Un
  i <- 1
  for(n in c(0:3)) {
    y.temp <- y.temp + ( Un[i]*(1-c)^(n+m)  ) / (1 + n*(1-E)/m)
    i <- i+1 #Update the counter variable
  } #next n
  
  # Next to complete the equation we multiply the summed component by the
  #   exploitation rate (E)
  
  # Lifetime yield from the cohort is ...
  y <- E*y.temp
  
  # To calculate yield-per-recruit we multiply lifetime yield by the mortality
  #   accrued between t0 and tr
  ypr <- y*exp(M*(tr-t0))*Winf
  
  # Return yield-per-recruit
  return(ypr)
}

# Lets test our function
BH_ypr(Fmort=0.1, M=0.1, tr=2,
                   c=0.5, # Ratio of Lc/Linf
                   Winf=2, k=0.2, t0=0)

BH_ypr(Fmort=0.5, M=0.1, tr=2,
       c=0.5, # Ratio of Lc/Linf
       Winf=2, k=0.2, t0=0)

# Now that we have this function yield-per-recruit in hand, lets 
#   simulate YPR across a range of fishing mortality rates (F) and values for c

trial.c <- seq(from=0.25, to=1, by=0.01)
trial.c

# Trial fishing mortality rates
trial.F <- seq(from=0, to=1, by=0.01)
trial.F

# To plot a surface of YPR as a function of these two parameters we 
#  can use expand_grid to create all combinations
?expand_grid

trial.df <- expand_grid(trial.c, trial.F)
head(trial.df)

# For additional perspective, lets calculate the exploitation rates
#  E=F/Z, that correspond to each of the F values 
#   assuming a natural mortality rate of M=0.1
M <- 0.1

# F=EZ
trial.df$trial.E<- trial.df$trial.F / ( trial.df$trial.F + M)

# OK, now that we have our trial values in the trial.df data frame,
head(trial.df)
tail(trial.df)

# We can loop through and calculate YPR for each combination of trial.c and trial.E
# First we will create an empty vector to hold the calculated YPR values
#   Of length equal to the length of our data frame of trial value combinations
n.trial.df <- nrow(trial.df)

YPR <- vector(length=n.trial.df)

# Now, we can loop through and calculate YPR with our function at each level of c and E
for(i in 1:n.trial.df) {
  YPR[i] <- BH_ypr(Fmort=trial.df$trial.F[i], M=M, tr=2,
                   c=trial.df$trial.c[i], # Ratio of Lc/Linf
                   Winf=2, k=0.2, t0=0)
}
YPR

# Now that we have calculated YPR for each of our trial values of F and c
#   lets add this variable to our data frame
trial.df$YPR <- YPR

# Yay! We made it, now the fun part. We will plot YPR as a function of c and F

# Plot the resulting surface
ggplot(data=trial.df, aes(x=trial.F, y=trial.c, z=YPR, fill=YPR)) +
  geom_raster( interpolate=TRUE, show.legend=TRUE) +
  scale_fill_viridis_c() +
  xlab("Instantaneous Fishing Mortality Rate (F)") +
  ylab("Ratio of Minimum Size Limit to Asymptotic Length (c = Lc/Linf)")

# Alternative plot:
ggplot(data=trial.df, aes(x=trial.F, y=trial.c, z=YPR)) +
  geom_contour_filled(show.legend=TRUE) +
  xlab("Instantaneous Fishing Mortality Rate (F)") +
  ylab("Ratio of Minimum Size Limit to Asymptotic Length (c = Lc/Linf)")

# Plot as points
ggplot(data=trial.df, aes(x=trial.F, y=trial.c, size=YPR, color=YPR)) +
  geom_point() +
  scale_color_viridis_c()


# If you want to see a cool alternative way to plot using the rayshader package
# install.packages("devtools")
# devtools::install_github("tylermorganwall/rayshader")
require(rayshader)

ggplt <- trial.df %>% 
  ggplot() +
  geom_tile(aes(x = trial.F, y = trial.c, fill = YPR)) +
  geom_contour(aes(x = trial.F, y = trial.c, z = YPR), color = "black") +
  scale_x_continuous("Fishing Mortality Rate (F)", expand = c(0, 0)) +
  scale_y_continuous("Size Limit Ratio (c= Lc/Linf)", expand = c(0, 0)) +
  scale_fill_gradientn("YPR", colours = terrain.colors(10)) +
  coord_fixed()

ggplt

# Now feed this into the plot_gg() function from rayshader
# WARNING: THIS MAY TAKE A WHILE TO RENDER! - JUST SKIP IF IT DOESN'T WORK!!!
plot_gg(ggplt, multicore = TRUE, raytrace = TRUE, width = 7, height = 4,
        scale = 300, windowsize = c(1400, 866), zoom = 0.6, phi = 30, theta = 30)

# Part D: Leslie Matrix Models =========================================
# Leslie (1945) showed that a linear discrete population model can be structured as a matrix equation,
#   where discrete time references (likely years) and discrete stages or ages are specified. 
# To do so we define a projection matrix whose elements describe the relationship between model stages.

# Numbers of individuals at age (or stage) are defined as a vector, and we can define our
#   projection matrix as we see fit.

# But first lets familiarize ourselves with matrix algebra operations in R.

# The %*% is how we specify a matrix multiplication operation. 

# Lets start with a simple example with vector of length 2, and a 2x2 projection matrix:

# Numbers at age (for example)
vect <- c(10,2)

# Projection matrix
mat <- matrix(data=c(5,10,
                     0.5,0), byrow=TRUE,
              nrow=2, ncol=2)

mat

# Here we can see what happens when we multiply our initial vector of number at age (vect)
#   by our projection matrix (mat)
vect.2 <- mat %*% vect
vect.2

# We can update our numbers at age, by multiply our new numbers at age vector (vect.2) again by our
#   projection matrix.
vect.3 <- mat %*% vect.2
vect.3

# NOTE, THE ORDER OF OPERATIONS MATTERS !!!!
mat %*% vect # Correct

# Is not the same as...
vect %*% mat # Incorrect

# OK, now that we have a handle on matrix multiplication within R, we can explore a more interesting
#   scenario.

# Here, we can define our projection matrix following a standard form,
#   where the elements in the first row correspond to the product of age-specific fecundities
#   and survival from birth to age-1.

# The off-diagonal elements that are not zero in our matrix represent our survival rates for a given age
#   survival for age-1 = 0.652, S for age-2 is 0.363, ect.

mat <- matrix(data=c(0,1.41,1.98,2.58,
                     0.652, 0, 0, 0,
                     0, 0.363, 0, 0,
                     0, 0, 0.128, 0),
                nrow=4, ncol=4, byrow=TRUE)
mat

# Now to project our numbers at age forward in time we need to define our initial numbers at age:
N0 <- c(1,1,1,1)

# Given this, our numbers at age after one year (first time step) will be
mat %*% N0

# Our numbers at age in the 2nd time step (after 2 years) can be found by multiplying our initial
#   vector of numbers at age (N0) by our projection matrix twice
mat %*% mat %*% N0

# To project our numbers at age further into the future, this becomes a bit repetitive.

# However, we can use the matrix.power() function in the "matrixcalc" package to
#   pre-multiply our projection matrix


mat %*% mat

matrix.power(mat,2)

# The following code will project our initial numbers at age (N0) two years forward in time...
# install.packages("matrixcalc")
library(matrixcalc)
matrix.power(mat,2) %*% N0

# This is great, but if we want to record our numbers at age in each time step,
#   we will need to loop through subsequent years and save the results to an object.

# Here we will project our age-structured population dynamics for 25 years, saving expected
#   numbers at age in each time step in a new matrix (N)

# 
n.years <- 25

# First, lets define a matrix to hold our numbers at age in each time step

N <- matrix(nrow=length(N0), ncol=n.years,
               dimnames=list(paste0("Age-",c(1:4)),
                             c(1:n.years)))

N

# Loop over forward projections
for(t in 0:(n.years-1)) {
  if(t==0) {
    N[,t+1] <- mat %*% N0 
  }else {
    N[,t+1] <- mat %*% N[,t]
  }
} #next t

N

# Now lets plot our resulting numbers at age

# For simplicity we will melt our matrix into a data frame for easy plotting with ggplot
# install.packages("reshape2")
library(reshape2)
N.df <- melt(N)
head(N.df)

# Lets give the data frame some useful names
names(N.df) <- c("Age","Year","N")

head(N.df)

# Now lets plot our expected numbers at age across time
g <- ggplot(N.df, aes(x=Year, y=N, fill=Age)) +
       theme_linedraw() +
       scale_fill_colorblind() +
       geom_area(alpha=0.75)
g

# If we zoom in a bit, we see that there is a period of years where the ratio
#   of ages is not constant.
g.sub <- g + scale_x_continuous(limits=c(1,10))
g.sub

# After this time period our population converges to what is known as the Stable Age Distribution (SAD)

# To explore different population dynamics, lets reduce our expected age-specific survival rates

# Alternative with population decline ============
mat.2 <- matrix(data=c(0,1.41,1.98,2.58,
                     0.4, 0, 0, 0,
                     0, 0.2, 0, 0,
                     0, 0, 0.1, 0),
              nrow=4, ncol=4, byrow=TRUE)
mat.2

# Again we will define a matrix to hold our expected numbers at age
n.years <- 25

N.2 <- matrix(nrow=length(N0), ncol=n.years,
            dimnames=list(paste0("Age-",c(1:4)),
                          c(1:n.years)))

# And our initial vector of numbers at age
N0.2 <- c(1e3,1e3,1e3,1e3)

for(t in 0:(n.years-1)) {
  if(t==0) {
    N.2[,t+1] <- mat.2 %*% N0.2
  }else {
    N.2[,t+1] <- mat.2 %*% N.2[,t]
  }
} #next t

N.2

# Now we can plot out our abundance at age over time
N.df.2 <- melt(N.2)
names(N.df.2) <- c("Age","Year","N")

N.df.2

g <- ggplot(N.df.2, aes(x=Year, y=N, fill=Age)) +
  theme_linedraw() +
  scale_fill_colorblind() +
  geom_area(alpha=0.75)
g

g.sub <- g + scale_x_continuous(limits=c(1,10))
g.sub

# EXPLORE THE IMPACT OF ASSUMING ALTERNATIVE SURVIAL RATES AND FECUNDITY-AGE0 SURVIVAL RATES...

# Part E: Pacific Salmon Matrix  =========================================

mat.salmon <- matrix(data=c(0, 0, 0, 50,
                            0.2, 0, 0, 0,
                            0, 0.3, 0, 0,
                            0, 0, 0.4, 0),
                     nrow=4, ncol=4, byrow=TRUE)

# Lets see the matrix
mat.salmon

N0.salmon <- c(1e4,1e3,1e2,1e2)

# Looping example
n.years <- 25

N.salmon <- matrix(nrow=length(N0.salmon), ncol=n.years,
                   dimnames=list(paste0("Ocean-",c(1:4)),
                             c(1:n.years)))

# As we saw before, we can either loop through years to calculate numbers at age
# with vector Nt+1 relative to Nt.

# Alternatively we can calculate Nt=2 as M*M*N0.

# We will use the max.power() function in the matrixcalc package
#   to conduct our matrix exponential
require(matrixcalc)
matrix.power(mat.salmon,2)

for(t in 1:n.years) {
  N.salmon[,t] <- matrix.power(mat.salmon,t) %*% N0.salmon
} #next t

N.salmon

# Convert simulated abundances to list for plotting
N.salmon.df <- melt(N.salmon)
names(N.salmon.df) <- c("Age","Year","N")

# Plot numbers at age across time
g <- ggplot(N.salmon.df, aes(x=Year, y=N, fill=Age)) +
  theme_linedraw() +
  scale_fill_colorblind() +
  geom_area(alpha=0.75) +
  ggtitle("Simplified Salmon Life History")
g

# For clarity, lets plot individual years
g <- ggplot(N.salmon.df, aes(x=Year, y=N, fill=Age)) +
  theme_linedraw() +
  scale_fill_colorblind() +
  geom_col(alpha=0.75) +
  ggtitle("Simplified Salmon Life History")
g

# WHAT TYPE OF DYNAMICS RESULT?
# WHAT IS THE INFLUENCE OF CHANGING FECUNDITY, OR EARLY OCEAN SURVIVAL?


# Part F: Stage-structured Matrix Models =========================================

# Leslie matrix models are not only useful in approximating age-structured population
#  dynamics, but are applicable to any finite time points in the species' life history.

# Here we can explore application of Leslie matrix model approach to life history stages.

# We will follow the example from Brewster-Giesz and Miller (2000) who
#   developed a life-stage model for sandbar shark.

# First, lets define our life stages
stages <- c("Neonates","Juveniles","Subadults","Adults","Resting")
n.stages <- length(stages)

# Next, we will define quantities within our projection matrix

# Probability of surviving the Neonate stage
gN <- 0.1

# Probability of surviving the Juvenile stage
gJ <- 0.3

# Probability of surviving the Subadult stage
gSA <- 0.5

# Probability of surviving the Adult stage
gA <- 0.8

# Fecundity for adults
FA <- 30

# Now, lets define some probabilities for remaining in the same stage
p.r.J <- 0.2 # Juveniles
p.r.SA <- 0.3 # Subadults

# Finally, the probability of moving from the Resting stage, back to the actively breeding
#   Adult stage
gR <- 0.5

# Now we will define our projection matrix
mat.shark <- matrix(data=c(0, 0, 0, gA*FA, 0,
                           gN, p.r.J, 0, 0, 0,
                           0, gJ, p.r.SA, 0, 0,
                           0, 0, gSA, 0, gR,
                           0, 0, 0, gA, 0),
                    nrow=n.stages, ncol=n.stages,
                    dimnames=list(stages,stages),
                    byrow=TRUE)

mat.shark

# Next, lets project this matrix forward across 25 years

N0.shark <- rep(1e3,n.stages)
N0.shark

n.years <- 25

N.shark <- matrix(nrow=n.stages, ncol=n.years,
                   dimnames=list(stages,
                                 c(1:n.years)))

# Leverage matrix multiplication to project abundance in each stage across years
require(matrixcalc)


for(t in 1:n.years) {
  N.shark[,t] <- matrix.power(mat.shark,t) %*% N0.shark
} #next t

N.shark

# Convert simulated abundances to list for plotting
N.shark.df <- melt(N.shark)
names(N.shark.df) <- c("Stage","Year","N")

# Plot number in each stage across time
g <- ggplot(N.shark.df, aes(x=Year, y=N, fill=Stage)) +
  theme_linedraw() +
  scale_fill_colorblind() +
  geom_area(alpha=0.75) +
  ggtitle("Shark Life History")
g

# For clarity, lets plot individual years
g <- ggplot(N.shark.df, aes(x=Year, y=N, fill=Stage)) +
  theme_linedraw() +
  scale_fill_colorblind() +
  geom_col(alpha=0.75) +
  ggtitle("Shark Life History")
g

# PLEASE EXPLORE HOW CHANGING YOUR THE VALUES OF THESE MODEL PARAMETERS ALTERS YOUR EXPECTATION.

















