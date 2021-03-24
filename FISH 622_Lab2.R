#==================================================================================================
#Project Name: FISH 622 - QUANTITATIVE FISH POPULATION DYNAMICS: Lab #2
#Creator: Curry James Cunningham, College of Fisheries and Ocean Sciences, UAF
#Date: January 29, 2021
#
#Purpose: Halibut CPUE Standardization
#
# =========================================================================================
require(tidyverse)
require(visreg)


# First, set your working directory!!!!!
setwd("/Users/curryc2/Documents/Students/2021/622 - Quant Fish Pop Dynamics/Content/Week 2/lab")


# Load CPUE Data ==============================================
cpue.dat <- read.csv(file="Halibut CPUE.csv", header=TRUE)

# Lets check it was read in correctly..
head(cpue.dat)
str(cpue.dat)


# Define Factors ==============================================
# Given that our Area and Gear predictors were read in as character strings, and
#   our year reference was read in as an integer, we will need to convert all of these to factors.

cpue.dat$fArea <- factor(cpue.dat$Area)
cpue.dat$fGear <- factor(cpue.dat$Gear)
cpue.dat$fYear <- factor(cpue.dat$Year)

# Lets check this worked
str(cpue.dat)

# Define reference level for factors ==========================
# Next we need to define the reference level for our factors, this will be:
#   Year: 1980
#   Area: SE-AK-Inside (US)
#   Gear: Fixed

# We can define the reference level for factors in R using the relevel() function
?relevel

cpue.dat$fYear <- relevel(cpue.dat$fYear, ref="1980")
cpue.dat$fArea <- relevel(cpue.dat$fArea, ref="SE-AK-Inside (US)")
cpue.dat$fGear <- relevel(cpue.dat$fGear, ref="Fixed")

# Fit GLM =====================================================
# Next we can fit our GLM

fit.cpue <- glm(ln_U ~ fGear + fYear + fArea, data=cpue.dat, family=gaussian(link="identity"))

# Lets see the summary of our estimated coefficients
summary(fit.cpue)

# How would you interpret these coefficients?

# Next lets visualize the effects
visreg(fit.cpue, type="contrast")


# Plot Resulting Fit =========================================
# To plot the resulting fit, lets add a new column to our data frame with the predicted ln(CPUE) i.e. pred_ln_U

# Predicted values
pred_ln_U <- predict(fit.cpue)

# Attach to data frame
cpue.dat <- data.frame(cpue.dat, pred_ln_U)

str(cpue.dat)

# Extract Coefficients
coef(fit.cpue)

# Now that we have both observed and predicted values, lets plot the fit

ggplot(data=cpue.dat, aes(x=Year, y=ln_U, color=fGear)) +
    theme_linedraw() +
    geom_point() +
    facet_wrap(~fArea) +
    geom_line(aes(y=pred_ln_U))

# Write out Predicted CPUE ====================================
write.csv(cpue.dat, "Updated Halibut CPUE.csv")
