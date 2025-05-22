install.packages(graphicalExtremes)
library("graphicalExtremes")
library(mev) # generating of multivariate extreme-value data
library(gtools) # calculating permutations
library(tailDepFun) # defining a grid
library(ggplot2) # plotting
library(tidyverse)
library(parallel)
main_application <- function( data , lambda_grid , grid  , num_col , start = NULL,
                              type = "SSR_row_log" , k, p, num_class = num_class, cl) {
  d <- ncol(data)  
  q <- nrow(grid)
  l <- d * num_col
  #checkN <- N %% num_class
  #if (checkN > 0) {
  #  N <- N - checkN
  #  cat("Warning: sample size N has been reduced to", N, "\n")
  #}
  
  # Data generation and preparation
  
  R <- apply(data, 2, rank)
  #Calculate the empirical stdf estimates on the grid points for the cross validation part
  w <- W_calculus(k = k, num_class = num_class, X = data, grid = grid, q = q)
  #Calculate the empitical stdf estimates on the grid points for the estimation part
  w_total <- sapply(1:q, function(m) stdfEmp(R, k, grid[m, ]))
  
  if (is.null(start)) {
    start <- starting_point(data, num_col)
    start <- c(t(start))
  }
  
  # Cross-validation wrapper
  wrapper_cross_validation <- function(lambda) {
    cross_validation_application(d, grid, lambda , num_col ,  start , type = type , p , w , num_class=num_class)
  }     
  
  # Parallelize cross-validation
  scores <- unlist(parLapply(cl, lambda_grid, wrapper_cross_validation))
  #print(scores)
  # Determine optimal lambda
  index <- which(scores == min(scores))
  lambda_optim <- lambda_grid[index]
  
  
  # Wrapper for param_estim
  wrapper_param_estim <- function(lambda) {
    param_estim_application(d , grid , lambda , num_col , start , type = type , p ,  w_total )
  }  
  
  # Execute param_estim in parallel
  Estimation <- parLapply(cl, as.list(lambda_optim), wrapper_param_estim)[[1]]
  # Compile results
  result <- list(lambda_optim = lambda_optim, Estimation = Estimation )
  
  return(result)
}



main_oversteps_application <- function( data , lambda_grid, grid,  start , 
                                        type = "SSR_row_log", k, p, num_class , cl ){
  num_col <- 1
  
  Estimation <- main_application( data , lambda_grid , grid  , num_col , start = start,
                                  type = "SSR_row_log" , k, p, num_class = num_class, cl)
  test <- TRUE
  while(test == TRUE){
    #print(num_col)
    old_Estimation <- Estimation
    old_matrix <- Estimation$Estimation$matrix
    sig_old <- apply(old_matrix, 2,  function(vec)   which(vec>0), simplify = FALSE )
    print(paste("The old estimation with", num_col, " number of columns is "))
    print(sig_old)
    print(old_Estimation)
    num_col <- num_col + 1 
    Estimation <- main_application( data , lambda_grid , grid  , num_col , start = start,
                                    type = "SSR_row_log" , k, p, num_class = num_class, cl)
    estim_matrix <- Estimation$Estimation$matrix
    sig_new  <- apply( estim_matrix, 2,  function(vec)   which(vec>0), simplify = FALSE )
    # print(sig_new)
    #print(paste("The new estimation with", num_col, " number of columns is ", sig_new))
    #if( has_duplicates(sig_new)  ) 
    #{ print(test)
    #test <- FALSE} 
    test <- all( apply( estim_matrix , 2 , function(vec) all(vec==0) ) == FALSE ) 
    print(test)
    
  }
  
  return(old_Estimation)
  
}










X <- danube$data_clustered[ , c(7 , 18 , 24 , 27 , 30) ]
X <- as.matrix(X)
rownames(X) <- NULL
colnames(X) <- NULL
d <- ncol(X)


# X: your n×d data matrix (no missing values)
n <- nrow(X)




#centers <-  c(0.04)
#generate_vector <- function(center, length = 20, range_factor = 0.4) {
#  step <- range_factor * center / (length / 2)  # Define step size proportionally
#  seq(from = center - (length / 2) * step + step / 2, 
#      to = center + (length / 2) * step - step / 2, 
#      length.out = length)
#}




####To generate the grid 

set.seed(123)  # for reproducibility

n_vectors <- 400
dim_vector <- 31
sparsity <- 0.8  # 80% zeros in each vector on average

# Function to generate one sparse vector
generate_sparse_vector <- function(d, sparsity) {
  x <- rep(0, d)
  num_nonzero <- ceiling((1 - sparsity) * d)
  indices <- sample(1:d, num_nonzero)
  x[indices] <- runif(num_nonzero, 0, 1)
  return(x)
}

# Generate matrix of 400 sparse vectors (each row is one vector)
#Grid_points_log <- t(replicate(n_vectors, generate_sparse_vector(dim_vector, sparsity)))








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

saveRDS(res0.4 , file = "C:/Users/mourahib/Desktop/github/Estimation-mixture-model/Results/application/five stations/res0.4.rds")

