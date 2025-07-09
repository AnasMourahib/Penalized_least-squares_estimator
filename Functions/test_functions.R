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
n <- 100
list <- list() 
count <- 0 
for(i in 1: n){          
  file_path <- paste0("C:/Users/mourahib/Desktop/github/Penalized_least-squares_estimator/Results/DA_mix_HR/extr_dir_ident/wrong_model/kn0.08/result_",i,".rds")
  if(file.exists(file_path)){
    count <- count + 1
    result <- readRDS(file_path)
    list[[count]] <- result 
  }
  else(print(i))
}
saveRDS(list , file= "C:/Users/mourahib/Desktop/github/Penalized_least-squares_estimator/Results/DA_mix_HR/extr_dir_ident/wrong_model/kn0.08.rds" )




###Effect of the tail fraction on teh score EDS      
n <- 100   #Replace p=0.2 by p=0.4 or p=0.6  4''
list <- readRDS(file= paste0("C:/Users/mourahib/Desktop/github/Penalized_least-squares_estimator/Results/DA_mix_log/effect_tail_fraction/N3000/pen0.2.rds"))
list[[1]] <- list[[2]]
#Calculus of EDS score. 
EDscore_p0.2 <- rep(0 , 10)
for( j in 1:10  ){
  list_matrices <- list()  
  for(i in 1: n){
    list_matrices <- list[[i]][[j]]$Estimation$matrix
    scorei <- EDS( list_matrices , A)
    EDscore_p0.2[j] <- EDscore_p0.2[j] + scorei/n
  }
}


#############plot EDS
colors <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00")  # Orange, Blue, Green

# X-axis values
x_values <- seq(0.02 , 0.2 , by=0.02) 

# Set up an empty plot with appropriate axis limits
par(mgp = c(2, 0.7, 0))  # Adjust label position
plot(x_values, EDscore_p0.2 , type = "l", col = colors[1], lty = 1, lwd = 2.5, 
     ylim =  c(0.5 , 1), 
     xlab = "",
      ylab = "" , #  "ED-S", 
     cex.lab = 1.5, cex.axis = 1.3, 
     xaxt = "n")

# Add custom x-axis ticks and labels
axis(1, at = x_values, labels = x_values, cex.axis = 1.3)




# Add lines for the other vectors
lines(x_values, EDscore_p0.4 , col = colors[2], lty = 2, lwd = 2.5)
lines(x_values, EDscore_p0.6 , col = colors[3], lty = 3, lwd = 2.5)
lines(x_values, EDscore_p1 , col = colors[4], lty = 3, lwd = 2.5)

# Add circles around the data points
points(x_values, EDscore_p0.2, col = colors[1], pch = 16, cex = 1.5)  # Solid circles
points(x_values, EDscore_p0.4 , col = colors[2], pch = 16, cex = 1.5)
points(x_values, EDscore_p0.6 , col = colors[3], pch = 16, cex = 1.5)
points(x_values, EDscore_p1 , col = colors[4], pch = 16, cex = 1.5)


###Effect of the tail fraction on the score SMSE 
n <- 100   #Replace p=0.2 by p=0.4 or p=0.6  4
list <- readRDS(file= paste0("C:/Users/mourahib/Desktop/github/Penalized_least-squares_estimator/Results/DA_mix_log/effect_tail_fraction/N3000/pen1.rds"))
list[[1]] <- list[[2]]
#Calculus of EDS score. 
list_res <- list()
SMSE_p1 <- rep(0 , 10)
for( j in 1:10  ){
  list_res <- list()
  for(i in 1 : n){
    list_res[[i]] <- list[[i]][[j]]$Estimation
  }
  SMSE_p1[j] <- SMSE_log_2(list_res , A , dep = 0.25)
}



#############plot EDS
colors <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00")  # Orange, Blue, Green

# X-axis values
x_values <- seq(0.02 , 0.2 , by=0.02) 

# Set up an empty plot with appropriate axis limits
par(mgp = c(2, 0.7, 0))  # Adjust label position
plot(x_values, SMSE_p0.2 , type = "l", col = colors[1], lty = 1, lwd = 2.5, 
     ylim =  c(0 , 4), 
     xlab = "tail fraction k/n",
     ylab = "" , #  "ED-S", 
     cex.lab = 1.5, cex.axis = 1.3, 
     xaxt = "n")

# Add custom x-axis ticks and labels
axis(1, at = x_values, labels = x_values, cex.axis = 1.3)




# Add lines for the other vectors
lines(x_values, SMSE_p0.4 , col = colors[2], lty = 2, lwd = 2.5)
lines(x_values, SMSE_p0.6 , col = colors[3], lty = 3, lwd = 2.5)
lines(x_values, SMSE_p1 , col = colors[4], lty = 3, lwd = 2.5)

# Add circles around the data points
points(x_values, SMSE_p0.2, col = colors[1], pch = 16, cex = 1.5)  # Solid circles
points(x_values, SMSE_p0.4 , col = colors[2], pch = 16, cex = 1.5)
points(x_values, SMSE_p0.6 , col = colors[3], pch = 16, cex = 1.5)
points(x_values, SMSE_p1 , col = colors[4], pch = 16, cex = 1.5)





#######Extreme directions identification 



############################Ectreme direction identification 

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




list <- readRDS(file= paste0("C:/Users/mourahib/Desktop/github/Penalized_least-squares_estimator/Results/DA_mix_HR/extr_dir_ident/wrong_model/kn0.08.rds"))
n <- 95

EDscore_kn0.08 <- rep(0  , 10)
for( j in 1:10  ){
  list_matrices <- list()
  for(i in 1: n){
    if(is.null(list[[i]]) == FALSE){
      list_matrices[[i]] <- list[[i]][[j]]$Estimation$pls_matrix
      scorei <- EDS( list_matrices[[i]] , A)
      if(scorei>1){print(3)}
      EDscore_kn0.08[j] <- EDscore_kn0.08[j] + scorei/n 
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
plot(x_values, EDscore_kn0.04 , type = "l", col = colors[1], lty = 1, lwd = 2.5, 
     ylim =  c(0.6 , 1), 
     xlab =  "n" , 
     ylab = "", 
     cex.lab = 1.5, cex.axis = 1.3, 
     xaxt = "n")

# Add custom x-axis ticks and labels
axis(1, at = x_values, labels = x_values, cex.axis = 1.3)




# Add lines for the other vectors
lines(x_values, EDscore_kn0.06 , col = colors[2], lty = 2, lwd = 2.5)
lines(x_values, EDscore_kn0.08, col = colors[3], lty = 3, lwd = 2.5)

# Add circles around the data points
points(x_values, EDscore_kn0.04, col = colors[1], pch = 16, cex = 1.5)  # Solid circles
points(x_values, EDscore_kn0.06 , col = colors[2], pch = 16, cex = 1.5)
points(x_values, EDscore_kn0.08 , col = colors[3], pch = 16, cex = 1.5)




list <- readRDS(file= paste0("C:/Users/mourahib/Desktop/github/Penalized_least-squares_estimator/Results/DA_mix_HR/effect_tail_fraction/N3000/pen0.2.rds"))
j <- 3
n <- 50
list_mat <- list()
for(i in 1 : n){
  list_mat[[i]] <- list[[i]][[j]]$Estimation$pls_dep
}

list_vec <- lapply(list_mat , matrix_vector )

my_list_trimmed <- lapply(list_vec <- lapply(list_mat , matrix_vector )
, function(x) x[-3])

# Step 2: Convert to a matrix (rows = observations, columns = remaining elements)
mat <- do.call(rbind, my_list_trimmed)

# Step 3: Define Greek-style labels using expression()
labels <- expression(Gamma[12], Gamma[13], Gamma[14], Gamma[24], Gamma[34])

# Step 4: Plot
bp <- boxplot(mat,
        main = "Sample size N = 3000",
        xlab = "Elements",
        ylab = "Values",
        names = labels, 
        ylim= c(0,5))
for (i in 1:length(bp$names)) {
  # You can tweak the width (here: 0.3 on each side of the box center)
  segments(x0 = i - 0.4, x1 = i + 0.4, y0 = 1, y1 = 1, col = "red", lty = 1, lwd = 2)
}


###########Effect of the penaliation exponent 

#Calculus of SMSE2 score
list <- readRDS(file= paste0("C:/Users/mourahib/Desktop/github/Penalized_least-squares_estimator/Results/DA_mix_log/effect_alpha/N3000/kn0.6.rds"))
list_estim <- list()
SMSE_kn0.6 <- rep(0 , 10)
for(j in 1 : 10){
  for( i in 1: n){
    list_estim[[i]] <- list[[i]][[j]]$Estimation
  }
  SMSE_kn0.6[j] <- SMSE_log_2(list_estim , A , dep = alpha )
}

#####plot SMSE

colors <- c("#E69F00", "#56B4E9", "#009E73")  # Orange, Blue, Green

# X-axis values
x_values <- seq(0.1 , 1 , by=0.1) 

# Set up an empty plot with appropriate axis limits
par(mgp = c(2, 0.7, 0))  # Adjust label position
plot(x_values, SMSE_kn0.4 , type = "l", col = colors[1], lty = 1, lwd = 2.5, 
     ylim =  c(0 , 1.2),    
     xlab = expression(p)
     , ylab = "" #"SMSE"
     , 
     cex.lab = 1.5, cex.axis = 1.3, 
     xaxt = "n")

# Add custom x-axis ticks and labels
axis(1, at = x_values, labels = x_values, cex.axis = 1.3)




# Add lines for the other vectors
lines(x_values, SMSE_kn0.5 , col = colors[2], lty = 2, lwd = 2.5)
lines(x_values, SMSE_kn0.6 , col = colors[3], lty = 3, lwd = 2.5)

# Add circles around the data points
points(x_values, SMSE_kn0.4 , col = colors[1], pch = 16, cex = 1.5)  # Solid circles
points(x_values, SMSE_kn0.5 , col = colors[2], pch = 16, cex = 1.5)
points(x_values, SMSE_kn0.6 , col = colors[3], pch = 16, cex = 1.5)








############################################################
##  3‑variate sample with mixed signs of dependence
##  X1–X3  positive   |   X1–X2, X2–X3  negative
############################################################

# install.packages("copula")   # run once if not installed
library(copula)

set.seed(42)                  # reproducibility
n <- 5000                     # sample size (rows of output)

## 1.  Build a correlation matrix with the required signs
rho12 <- -0.40   # ρ(X1,X2)  < 0
rho13 <-  0.60   # ρ(X1,X3)  > 0
rho23 <- -0.40   # ρ(X2,X3)  < 0

R <- matrix(c( 1 ,  rho12, rho13,
               rho12, 1   , rho23,
               rho13, rho23, 1   ), nrow = 3, byrow = TRUE)

# Check positive‑definiteness (should be TRUE)
eigen(R)$values > 0        # all TRUE?  good!

## 2.  Define the Gaussian copula with that correlation
gauss_cop <- normalCopula(param = c(rho12, rho13, rho23),
                          dim   = 3, dispstr = "un")

## 3.  Draw n samples of uniforms (U1,U2,U3) from the copula
u <- rCopula(n, gauss_cop)        # n × 3 matrix

## 4.  Optional: map to any margins you like
#     (here we keep them Uniform(0,1); change with qnorm(), qexp(), etc.)

## 5.  Quick verification of signs via sample Kendall’s tau
kendall <- cor(u, method = "kendall")
kendall
#          [,1]      [,2]      [,3]
# [1,]  1.00000 -0.272xxx  0.411xxx   (signs match design)
# [2,] -0.272xxx  1.00000 -0.271xxx
# [3,]  0.411xxx -0.271xxx  1.00000

## 6.  Your simulated data matrix
head(u)


lambda_grid <- seq(0.01, 1 , length.out = 50)
p <- 0.4
k <- nrow(sim_data)/10
num_class <- 5  
points_log <- c(0,1/4 , 1/3 , 1/2 , 3/4 ,1)
Grid_points_log <- selectGrid(cst = points_log, d = 3, nonzero  = c( 2,3) )

cl <- makeCluster(6)
# Export necessary functions and variables to cluster
clusterExport(cl, varlist = c(   "main_application", "shuffleCols" ,"param_estim_application", "construct_symmetric_matrix",  "normalize_group", "cross_validation_application", "p", "k","Grid_points_log" , "X"  , "num_class"
                                 , "lambda_grid" 
))
clusterEvalQ(cl, { dyn.load("main.dll") })

set.seed(123) 





res <- main_oversteps_application( data = u , lambda_grid , grid = Grid_points_log,  start = NULL , 
                                   type = "SSR_row_log", k, p, num_class , cl )

