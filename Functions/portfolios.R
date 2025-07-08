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


data <- read.csv("C:/Users/mourahib/Desktop/github/Penalized_least-squares_estimator/Data/10_Industry_Portfolios_Daily.csv",
                 skip = 9,
                 stringsAsFactors = FALSE)



# Keep only the 10 return columns (drop date and row number)
returns_matrix <- as.matrix(data[, 2:11])
colnames(returns_matrix) <- NULL
returns_matrix <- matrix(as.numeric(returns_matrix), nrow = nrow(returns_matrix), ncol = ncol(returns_matrix))
bad_rows <- which(!complete.cases(mydata))   # positions of rows with NA or NaN
returns_matrix <- returns_matrix[ -bad_rows ,]
log_returns <-  log((returns_matrix/100) + 1 )
data <- - log_returns

#We fit a GARCH(1,1) model to each margin as in Equation (5.4) and extract the residuals to remove time dependency 

eps_list <- lapply(1:ncol(data), function(i) {
  gfit <- fit_GARCH_11(data[, i])
  gfit$Z.t
})
eps_matrix <- do.call(cbind, eps_list)


data <- eps_matrix

d <- ncol(data)
###For p= 0.4, change lambda_grid to : seq(0.001, 0.02 ,  by = 0.002)
lambda_grid <- seq(0.035, 0.055 ,  by = 0.002)
p <- 0.6
k <- nrow(X)/20
num_class <- 5  
points_log <- c(0,1/3  , 2/3 ,1)
Grid_points_log <- selectGrid(cst = points_log, d = d, nonzero  = c(1, 2) )




cl <- makeCluster(6)
# Export necessary functions and variables to cluster
clusterExport(cl, varlist = c(   "main_application", "shuffleCols" ,"param_estim_application", "construct_symmetric_matrix",  "normalize_group", "cross_validation_application", "p", "k","Grid_points_log" , "X"  , "num_class"
                                 , "lambda_grid" 
))
clusterEvalQ(cl, { dyn.load("main.dll") })

set.seed(123) 



res0.6 <- main_oversteps_application( data = data , lambda_grid , grid = Grid_points_log,  start = NULL , 
                                      type = "SSR_row_log", k, p, num_class , cl )
saveRDS(res0.6, file = "C:/Users/mourahib/Desktop/github/Penalized_least-squares_estimator/Results/application/portfolios10res0.6.rds")

########### Extracting the results 


res <- readRDS(file = "C:/Users/mourahib/Desktop/github/Penalized_least-squares_estimator/Results/application/portfolios10res0.4.rds")
#mat is the estimated matrix and coeff the estimated dependence coefficient for the mixture logistic model 
mat <- res$Estimation$matrix
coeff <- res$Estimation$dep
#extreme directions is a list of the resulting extreme directions
extreme_directions <- list()
for(i in 1 : ncol(mat))
{
  extreme_directions[[i]] <- which(mat[,i] >0  )
}
lextreme_direction <- length(extreme_directions)
#We calcualte the weight of each extreme direction as in Equation~(5.2) from the paper  
f <- function(vec){
  sum( vec  )^(coeff)
}
mat_pow <- mat^(1/coeff)
pre_weights <- apply(mat_pow , 2 , f   )
mass_tot <-  sum( pre_weights ) 
weights <-  pre_weights / mass_tot
