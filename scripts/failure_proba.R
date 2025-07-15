library(parallel)
library(graphicalExtremes)
library(mev) # generating of multivariate extreme-value data
library(gtools) # calculating permutations
library(tailDepFun) # defining a grid
library(ggplot2) # plotting
library(tidyverse)
dyn.load("main.dll")
source("Functions/Simulation_mixture_model.R")
source("Functions/Starting_points.R")
source("Functions/Help_functions.R")
source("Functions/param_estim.R")
source("Functions/cross_validation.R")
source("Functions/main.R")
source("Functions/main_applications.R")






#####Extracting data


X <- danube$data_clustered[ , c(7 , 18 , 24 , 27 , 30) ]
X <- as.matrix(X)

d <- ncol(X)
n <- nrow(X)

Max_per_observ <- apply(X ,  1 , max  )
#q0.9 is the 0.9 quantile
q0.9 <- quantile(Max_per_observ , probs = 0.9)
#Now, we calculate the failure probability using the model fit

#Calculate the componentwise maxima each year # Assume X is your matrix
years <- rownames(X)  # extract year labels
X_split <- split(data.frame(X), years)  # split by year (turn matrix into data.frame first)

# Compute column-wise max for each year
yearly_max <- t(sapply(X_split, function(df) apply(df, 2, max)))
yearly_max <- as.matrix(yearly_max)
rownames(yearly_max) <- NULL
colnames(yearly_max) <- NULL



