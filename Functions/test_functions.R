library(parallel)
library(mev) # generating of multivariate extreme-value data
library(gtools) # calculating permutations
library(tailDepFun) # defining a grid
library(ggplot2) # plotting
library(tidyverse)
install.packages("profmem")  # Only needed the first time
library(profmem)
dyn.load("main.dll")
source("Functions/Simulation_mixture_model.R")
source("Functions/Starting_points.R")
source("Functions/Help_functions.R")
source("Functions/param_estim.R")
source("Functions/cross_validation.R")
source("Functions/main.R")





N <- 10^6
d <- 4 
r <- 4 
num_col <- 4 
points <- c(0 , 1/6,  1/8 ,  1/4 ,  1/3 , 1/2 , 2/3 , 3/4 , 1)
Grid_points_HR <- selectGrid(points, d = d, nonzero = 2  )
lambda <- 10^(-3)
lambda_sq <- 1/8
Sigma <- matrix(lambda_sq , nrow = d , ncol = d )
diag(Sigma) <- 0
X <- N_generate_Mix_hr( N = N , A , Sigma  )
R <- apply(X , 2 , rank)
start <- starting_point(X, num_col)
start <- c(t(start))
p <- 0.3
q <- nrow(Grid_points_HR)
k <- 0.1 
w_total <- sapply(1:q, function(m) stdfEmp(R, k * N, Grid_points_HR[m, ]))


#####Test for param_estim 

param_estim(d , r , A , grid = Grid_points_HR , lambda , num_col = NULL , start , type = "SSR_row_HR", p = p ,  w = w_total )

##Test for cross_validation
w <- W_calculus(k = k * N , num_class = 10 , X = X, grid = Grid_points_HR, q = q)
cross_validation(d , r, A , grid = Grid_points_HR, lambda = 0.001 , num_col = NULL, start , type = "SSR_row_HR", p = 0.3 , w , num_class=10)


#### Test for main function
lambda_grid <- lambda_grid[[1]]
N <- 10^6




lambda_grid <- 1e-06


f <- function(seed){
  cl <- makeCluster(7)
  # Export necessary functions and variables to cluster
  clusterExport(cl, varlist = c("shuffleCols" ,"param_estim", "construct_symmetric_matrix",  "normalize_group", "cross_validation", "p", "d", "r", "Grid_points_HR"  , "num_col" , "A", "permutations"))
  clusterEvalQ(cl, { dyn.load("main.dll") })  # Or load "main.so" depending on your system
 res <-  main(lambda_grid, N, seed, A, grid = Grid_points_HR, alpha = NULL, Sigma = Sigma ,  num_col = NULL, start = NULL,
       type = "SSR_row_HR", k = N * k, p = 0.1 , num_class = 10, cl)
 stopCluster(cl)
 return(res)
}

mem_used <- profmem(res <- f(seed = 123))


list2 <- list()
for(i in 1 : 1){
  print(i)
  list2[[i]] <- f(i)
}











