library(parallel)
library(graphicalExtremes)
library(mev) # generating of multivariate extreme-value data
library(gtools) # calculating permutations
library(tailDepFun) # defining a grid
library(ggplot2) # plotting
library(tidyverse)
library(vars)
library(qrmtools)
dyn.load("main.dll")
source("Functions/Simulation_mixture_model.R")
source("Functions/Starting_points.R")
source("Functions/Help_functions.R")
source("Functions/param_estim.R")
source("Functions/cross_validation.R")
source("Functions/main.R")
source("Functions/main_applications.R")

#########Ewtracting data


########Simulation
simu_frechet <- function(n){
  U <- runif(n , min = 0 , max = 1)
  X = -1/log(U)
  return(X)
}


empirical_cdf <- function(x) {
  # rank each x[i] by how many values are ≤ it,
  # then divide by n to get the empirical CDF at x[i]
  rank(x, ties.method = "max") / length(x)
}

uniform_Frechet_transform <- function(U){
  X = -1/log(U)
  return(X)
}

n = 10^4
set.seed(7)
N1 = rnorm(n , 0 , 1)
N2 =    rnorm(n , 0 , 1)


X1 = N1
X2 = X1^2 + N2 #apply( cbind(N1 , N2), 1 ,  max)

data = cbind(X1 , X2)

###########Tuning parameters
generate_vector <- function(center, length = 20, range_factor = 0.4) {
  step <- range_factor * center / (length / 2)  # Define step size proportionally
  seq(from = center - (length / 2) * step + step / 2,
      to = center + (length / 2) * step - step / 2,
      length.out = length)
}
lambda_grid <- seq(0.001 , 0.01 , by = 0.001)
p <- 0.4
k <- nrow(data)/20
num_class <- 10
points_log <- c(0,1/3  , 2/3 ,1)
d = 2
Grid_points_log <- selectGrid(cst = points_log, d = d, nonzero  = c( 2) )


############Fit



cl <- makeCluster(6)
# Export necessary functions and variables to cluster
clusterExport(cl, varlist = c(   "main_application", "shuffleCols" ,"param_estim_application", "construct_symmetric_matrix",  "normalize_group", "cross_validation_application", "p", "k","Grid_points_log" , "data"  , "num_class"
                                 , "lambda_grid"
))
clusterEvalQ(cl, { dyn.load("main.dll") })

set.seed(123)



res <- main_oversteps_application( data = data , lambda_grid , grid = Grid_points_log,  start = NULL ,
                                   type = "SSR_row_log", k, p, num_class , cl )





################New data set graph with directed edges 1->2, 1 ->3, 2->4 and 3 ->4
a12 = 1
a13 = 1
a24 = 1
a34 = 1

n = 10^5
set.seed(7)
N1 = rnorm(n , 0 , 1)
N2 = rnorm(n , 0 , 1)
N3 = rnorm(n , 0 , 1)
N4 = rnorm(n , 0 , 1)


X1 = N1
X2 = apply(cbind(a12 * X1  , N2), 1, max )
X3 = apply(cbind(a13 * X1  , N3), 1, max )
X4 = apply(cbind(a24 * X2 , a34 * X3 , N4), 1, max )


data = cbind(X1 , X2 , X3 , X4)

###########Tuning parameters
generate_vector <- function(center, length = 20, range_factor = 0.4) {
  step <- range_factor * center / (length / 2)  # Define step size proportionally
  seq(from = center - (length / 2) * step + step / 2,
      to = center + (length / 2) * step - step / 2,
      length.out = length)
}
lambda_grid <- generate_vector(0.01)
p <- 0.4
k <- nrow(data)/20
num_class <- 10
points_log <- c(0,1/3  , 2/3 ,1)
d = 4
Grid_points_log <- selectGrid(cst = points_log, d = d, nonzero  = c( 2, 3 , 4) )




############Fit



cl <- makeCluster(6)
# Export necessary functions and variables to cluster
clusterExport(cl, varlist = c(   "main_application", "shuffleCols" ,"param_estim_application", "construct_symmetric_matrix",  "normalize_group", "cross_validation_application", "p", "k","Grid_points_log" , "data"  , "num_class"
                                 , "lambda_grid"
))
clusterEvalQ(cl, { dyn.load("main.dll") })

set.seed(123)



res <- main_oversteps_application( data = data , lambda_grid , grid = Grid_points_log,  start = NULL ,
                                   type = "SSR_row_log", k, p, num_class , cl )



############# I. Do this in the case where casual structures appear on extremes
###########I.1 X1 -> X2 when X1 is large with max(X1 , X2) and X2 = X1^(1/2) + X2


########Simulation
simu_frechet <- function(n){
  U <- runif(n , min = 0 , max = 1)
  X = -1/log(U)
  return(X)
}


empirical_cdf <- function(x) {
  # rank each x[i] by how many values are ≤ it,
  # then divide by n to get the empirical CDF at x[i]
  rank(x, ties.method = "max") / length(x)
}

uniform_Frechet_transform <- function(U){
  X = -1/log(U)
  return(X)
}

n = 10^5
set.seed(7)
N1 = rnorm(n , 0 , 1)
N2 = rnorm(n , 0 , 1)
X1 = N1
X2 = N2
##########We define a threshold u, the q0.95 quantile of X1. When X1 is larger than this threshold, then there's a causal effect X1 -> X2
u = quantile(X1 , probs = 0.95 ) # for both "X2 = (X1)^(1/2) + X2" and max, this works with probs = 0.95
extreme_obs = which(X1 > u)
X2[extreme_obs] = (X1[extreme_obs])^(1/2) + X2[extreme_obs] ### apply(cbind(X2[extreme_obs] , X1[extreme_obs]) , 1 , max)  #(X1[extreme_obs])^(1/2) + X2[extreme_obs] #apply(cbind(X2[extreme_obs] , X1[extreme_obs]) , 1 , max)    #(X1[extreme_obs])^2 + X2[extreme_obs]   #apply(cbind(X2[extreme_obs] , X1[extreme_obs]) , 1 , max)


data = cbind(X1 , X2)


###########Tuning parameters
generate_vector <- function(center, length = 20, range_factor = 0.4) {
  step <- range_factor * center / (length / 2)  # Define step size proportionally
  seq(from = center - (length / 2) * step + step / 2,
      to = center + (length / 2) * step - step / 2,
      length.out = length)
}
lambda_grid <- generate_vector(0.003)  ### for the max, the lambda grid that works is centered around 0.03 with sample size n = 10^5. For "X2 = (X1)^(1/2) + X2", this works with a grid centered around 0.003
p <- 0.4
k <- nrow(data)/20 ### For both "max" and "X2 = (X1)^(1/2) + X2", this works with this works with k = nrow(data)/20
num_class <- 10
points_log <- c(0,1/3  , 1/4 , 1/2 ,  2/3 , 7/8 ,1)
d = 2
Grid_points_log <- selectGrid(cst = points_log, d = d, nonzero  = c( 2) )




############Fit



cl <- makeCluster(6)
# Export necessary functions and variables to cluster
clusterExport(cl, varlist = c(   "main_application", "shuffleCols" ,"param_estim_application", "construct_symmetric_matrix",  "normalize_group", "cross_validation_application", "p", "k","Grid_points_log" , "data"  , "num_class"
                                 , "lambda_grid"
))
clusterEvalQ(cl, { dyn.load("main.dll") })

set.seed(123)



res <- main_oversteps_application( data = data , lambda_grid , grid = Grid_points_log,  start = NULL ,
                                   type = "SSR_row_log", k, p, num_class , cl )





















###########I.2 X1 -> X2 when X1 is large with X2 = X1^(2) + X2 (Not working yet)


########Simulation
simu_frechet <- function(n){
  U <- runif(n , min = 0 , max = 1)
  X = -1/log(U)
  return(X)
}


empirical_cdf <- function(x) {
  # rank each x[i] by how many values are ≤ it,
  # then divide by n to get the empirical CDF at x[i]
  rank(x, ties.method = "max") / length(x)
}

uniform_Frechet_transform <- function(U){
  X = -1/log(U)
  return(X)
}

n = 10^5
set.seed(7)
N1 = rnorm(n , 0 , 1)
N2 = rnorm(n , 0 , 1)
X1 = N1
X2 = N2
##########We define a threshold u, the q0.95 quantile of X1. When X1 is larger than this threshold, then there's a causal effect X1 -> X2
u = quantile(X1 , probs = 0.95 ) # for both "X2 = (X1)^(1/2) + X2" and max, this works with probs = 0.95
extreme_obs = which(X1 > u)
X2[extreme_obs] = (X1[extreme_obs])^(2) + X2[extreme_obs] ### apply(cbind(X2[extreme_obs] , X1[extreme_obs]) , 1 , max)  #(X1[extreme_obs])^(1/2) + X2[extreme_obs] #apply(cbind(X2[extreme_obs] , X1[extreme_obs]) , 1 , max)    #(X1[extreme_obs])^2 + X2[extreme_obs]   #apply(cbind(X2[extreme_obs] , X1[extreme_obs]) , 1 , max)


data = cbind(X1 , X2)


###########Tuning parameters
generate_vector <- function(center, length = 20, range_factor = 0.4) {
  step <- range_factor * center / (length / 2)  # Define step size proportionally
  seq(from = center - (length / 2) * step + step / 2,
      to = center + (length / 2) * step - step / 2,
      length.out = length)
}
lambda_grid <- generate_vector(0.003)  ### for the max, the lambda grid that works is centered around 0.03 with sample size n = 10^5. For "X2 = (X1)^(1/2) + X2", this works with a grid centered around 0.003
p <- 0.4
k <- nrow(data)/20 ### For both "max" and "X2 = (X1)^(1/2) + X2", this works with this works with k = nrow(data)/20
num_class <- 10
points_log <- c(0,1/3  , 1/4 , 1/2 ,  2/3 , 7/8 ,1)
d = 2
Grid_points_log <- selectGrid(cst = points_log, d = d, nonzero  = c( 2) )




############Fit



cl <- makeCluster(6)
# Export necessary functions and variables to cluster
clusterExport(cl, varlist = c(   "main_application", "shuffleCols" ,"param_estim_application", "construct_symmetric_matrix",  "normalize_group", "cross_validation_application", "p", "k","Grid_points_log" , "data"  , "num_class"
                                 , "lambda_grid"
))
clusterEvalQ(cl, { dyn.load("main.dll") })

set.seed(123)



res <- main_oversteps_application( data = data , lambda_grid , grid = Grid_points_log,  start = NULL ,
                                   type = "SSR_row_log", k, p, num_class , cl )

