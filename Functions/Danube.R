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










X <- danube$data_clustered[ , c(7 , 18 , 24 , 27 , 30) ]
X <- as.matrix(X)
rownames(X) <- NULL
colnames(X) <- NULL
d <- ncol(X)



n <- nrow(X)















lambda_grid <- seq(0.01, 1 , length.out = 100)
p <- 0.4
k <- nrow(X)/10
num_class <- 5  
points_log <- c(0,1/3 , 1/2 , 2/3 ,1)
Grid_points_log <- selectGrid(cst = points_log, d = d, nonzero  = c(1, 2,3) )

cl <- makeCluster(6)
# Export necessary functions and variables to cluster
clusterExport(cl, varlist = c(   "main_application", "shuffleCols" ,"param_estim_application", "construct_symmetric_matrix",  "normalize_group", "cross_validation_application", "p", "k","Grid_points_log" , "X"  , "num_class"
                                 , "lambda_grid" 
))
clusterEvalQ(cl, { dyn.load("main.dll") })

set.seed(123) 
res0.4 <- main_oversteps_application( data = X , lambda_grid , grid = Grid_points_log,  start = NULL , 
                                      type = "SSR_row_log", k, p, num_class , cl )
