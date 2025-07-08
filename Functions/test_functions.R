library(parallel)
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




list2 <- list()
for(i in 1 : 1){
  print(i)
  list2[[i]] <- f(i)
}





#############Transform data 
n <- 50
list <- list() 
for(i in 1: n){         
  file_path <- paste0("C:/Users/mourahib/Desktop/github/Penalized_least-squares_estimator/Results/DA_mix_log/effect_tail_fraction/N10000/pen1/result_",i,".rds")
  if(file.exists(file_path)){
    result <- readRDS(file_path)
    list[[i]] <- result 
  }
  else(print(i))
}
saveRDS(list , file= "C:/Users/mourahib/Desktop/github/Penalized_least-squares_estimator/Results/DA_mix_log/effect_tail_fraction/N10000/pen0.4.rds" )




###Effect of the tail fraction on teh score EDS      
n <- 100  #Replace p=0.2 by p=0.4 or p=0.6  4''
list <- readRDS(file= paste0("C:/Users/mourahib/Desktop/github/Penalized_least-squares_estimator/Results/DA_mix_log/effect_tail_fraction/N3000/pen0.2.rds"))
#Calculus of EDS score. 
EDscore_N3 <- rep(0 , 10)
for( j in 1:10  ){
  list_matrices <- list()  
  for(i in 1: n){
    list_matrices <- list[[i]][[j]]$Estimation$matrix
    scorei <- EDS( list_matrices , A)
    EDscore_N3[j] <- EDscore_N3[j] + scorei/n
  }
}


#############plot EDS
colors <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00")  # Orange, Blue, Green

# X-axis values
x_values <- seq(0.02 , 0.2 , by=0.02) 

# Set up an empty plot with appropriate axis limits
par(mgp = c(2, 0.7, 0))  # Adjust label position
plot(x_values, EDscore_N1 , type = "l", col = colors[1], lty = 1, lwd = 2.5, 
     ylim =  c(0 , 0.6), 
     xlab =  "tail fraction k/n",
      ylab = "ED-S" , #  "ED-S", 
     cex.lab = 1.5, cex.axis = 1.3, 
     xaxt = "n")

# Add custom x-axis ticks and labels
axis(1, at = x_values, labels = x_values, cex.axis = 1.3)




# Add lines for the other vectors
lines(x_values, EDscore_N2 , col = colors[2], lty = 1, lwd = 2.5)
lines(x_values, EDscore_N3 , col = colors[3], lty = 1, lwd = 2.5)

# Add circles around the data points
points(x_values, EDscore_N1, col = colors[1], pch = 16, cex = 1.5)  # Solid circles
points(x_values, EDscore_N2 , col = colors[2], pch = 16, cex = 1.5)
points(x_values, EDscore_N3 , col = colors[3], pch = 16, cex = 1.5)



###Effect of the tail fraction on the score SMSE 
n <- 50   #Replace p=0.2 by p=0.4 or p=0.6  4''
list <- readRDS(file= paste0("C:/Users/mourahib/Desktop/github/Penalized_least-squares_estimator/Results/DA_mix_HR/effect_tail_fraction/N3000/pen0.4.rds"))
#Calculus of EDS score. 
SMSE_N3 <- rep(0 , 10)
for( j in 1:10  ){
  list_res <- list()
  for(i in 1 : n){
    list_res[[i]] <- list[[i]][[j]]
  }
  SMSE_N3[j] <-  SMSE_HR2 (list_res , A , ind = 3  ) # SMSE_log_2 (list_res , A , dep = 0.25)   #  
}



#############plot EDS
colors <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00")  # Orange, Blue, Green

# X-axis values
x_values <- seq(0.02 , 0.2 , by=0.02) 

# Set up an empty plot with appropriate axis limits
par(mgp = c(2, 0.7, 0))  # Adjust label position
plot(x_values, SMSE_N1 , type = "l", col = colors[1], lty = 1, lwd = 2.5, 
     ylim =  c(0 , 5), 
     xlab = "tail fraction k/n",
     ylab = "SMSE" , #  "ED-S", 
     cex.lab = 1.5, cex.axis = 1.3, 
     xaxt = "n")

# Add custom x-axis ticks and labels
axis(1, at = x_values, labels = x_values, cex.axis = 1.3)




# Add lines for the other vectors
lines(x_values, SMSE_N2 , col = colors[2], lty = 1, lwd = 2.5)
lines(x_values, SMSE_N3 , col = colors[3], lty = 1, lwd = 2.5)

# Add circles around the data points
points(x_values, SMSE_N1, col = colors[1], pch = 16, cex = 1.5)  # Solid circles
points(x_values, SMSE_N2 , col = colors[2], pch = 16, cex = 1.5)
points(x_values, SMSE_N3 , col = colors[3], pch = 16, cex = 1.5)






#######Extreme directions identification 



############################Extreme direction identification 

n <- 100
count <- 0 
list <- list() 
for(i in 1: n){         
  file_path <- paste0("C:/Users/mourahib/Desktop/github/Penalized_least-squares_estimator/Results/DA_mix_HR/extr_dir_ident/kn0.08/result_" , i, ".rds")
  if(file.exists(file_path)){
    count <- count + 1
    result <- readRDS(file_path)
    list[[count]] <- result 
  }
  else(print(i))
}
saveRDS(list , file= "C:/Users/mourahib/Desktop/github/Penalized_least-squares_estimator/Results/DA_mix_HR/extr_dir_ident/kn0.08.rds" )





n <- 100
for(i in 1 : n){
  file_path <-  paste0("C:/Users/mourahib/Desktop/github/Penalized_least-squares_estimator/Results/DA_mix_HR/extr_dir_ident/kn0.04/result_" , i, ".rds")
  if(file.exists(file_path) == FALSE){
    print(i)
  }
}




list <- readRDS(file= paste0("C:/Users/mourahib/Desktop/github/Estimation-mixture-model/Results/extreme_dir_iden/kn0.08.rds"))
n <- 100

EDscore_kn0.082 <- rep(0  , 10)
for( j in 1:10  ){
  list_matrices <- list()
  for(i in 1: n){
    if(is.null(list[[i]]) == FALSE){
      list_matrices[[i]] <- list[[i]][[j]]$Estimation$matrix
      scorei <- EDS( list_matrices[[i]] , A)
      if(scorei>1){print(3)}
      EDscore_kn0.082[j] <- EDscore_kn0.082[j] + scorei/n 
    }
    else(print(i))
  }
}

##########plot EDscore extr_iden
colors <- c("#E69F00", "#56B4E9", "#009E73")  # Orange, Blue, Green

# X-axis values
x_values <- seq(500 , 5000 , by = 500)

# Set up an empty plot with appropriate axis limits
par(mgp = c(2, 0.7, 0))  # Adjust label position
plot(x_values, EDscore_kn0.08 , type = "l", col = colors[1], lty = 1, lwd = 2.5, 
     ylim =  c(0 , 0.2), 
     xlab =  "sample size n" , 
     ylab = "ED-S", 
     cex.lab = 1.5, cex.axis = 1.3, 
     xaxt = "n")

# Add custom x-axis ticks and labels
axis(1, at = x_values, labels = x_values, cex.axis = 1.3)




# Add lines for the other vectors
lines(x_values, EDscore_kn0.082 , col = colors[2], lty = 1, lwd = 2.5)
#lines(x_values, EDscore_kn0.08, col = colors[3], lty = 1, lwd = 2.5)

# Add circles around the data points
points(x_values, EDscore_kn0.08, col = colors[1], pch = 16, cex = 1.5)  # Solid circles
points(x_values, EDscore_kn0.082 , col = colors[2], pch = 16, cex = 1.5)
points(x_values, EDscore_kn0.08 , col = colors[3], pch = 16, cex = 1.5)


