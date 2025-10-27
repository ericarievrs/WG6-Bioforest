
library(data.table)
library(ggplot2)
library(effects)
library(interactions)
library(tidyverse)


set.seed(123)

dat <- data.frame(
  Treat = factor(rep(c("Control", "Disturbed"), each = 60))
)

dat$CWM <- runif(120, 2, 5) 
dat$VPD <- runif(120, 1, 3)

# Define growth with different VPD sensitivity
dat$growth <- ifelse(dat$Treat == "Disturbed",
                     5 + 1.2*dat$CWM - 1.5*dat$VPD + rnorm(60, 0, 0.5),  # disturbed: strong VPD effect
                     5 + 1.2*dat$CWM - 0.5*dat$VPD + rnorm(60, 0, 0.5))  # control: weaker VPD effect

model <- lm(growth ~ VPD * CWM * Treat, data = dat)
summary(model)

plot(allEffects(model))

interact_plot(model, pred = "VPD", modx = "CWM", plot.points = F)

interact_plot(model,
              pred = VPD,
              modx = CWM,
              mod2 = Treat, 
              modx.values = "plus-minus",
              plot.points = FALSE,
              interval = FALSE) +
  theme_minimal()




