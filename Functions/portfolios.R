library(parallel)
library(graphicalExtremes)
library(mev) # generating of multivariate extreme-value data
library(gtools) # calculating permutations
library(tailDepFun) # defining a grid
library(ggplot2) # plotting
library(tidyverse)
library(vars)
dyn.load("main.dll")
source("Functions/Simulation_mixture_model.R")
source("Functions/Starting_points.R")
source("Functions/Help_functions.R")
source("Functions/param_estim.R")
source("Functions/cross_validation.R")
source("Functions/main.R")
source("Functions/main_applications.R")





tmp <- read.csv("C:/Users/mourahib/Desktop/github/Graph-estimation/data/30_Industry_Portfolios_Daily.CSV", header= T, nrows = 24643)
tmp$X <- as.Date(as.character(tmp$X), format = "%Y%m%d")
indx <- which(tmp$X == "1970-01-02") # Could pick another date but check for missing values
dates <- tmp$X[indx:length(tmp$X)]
data <- as.matrix(-tmp[indx:nrow(tmp),2:ncol(tmp)])
rownames(data) <- NULL
colnames(data) <- NULL




##Decluster data using an AR(p) model with p chosen using an AIC criterion 
## assume your data are in a matrix 'data' of dimension T × d
# assume your data are in a matrix 'data' of dimension T × d
data <- data[, c(1 , 2 ,3 , 4 , 5) ] ##we take only 8 portfolios



# Suppose your data are in a T×30 matrix called `data`

# 1) Use VARselect() to pick p by AIC
#    lag.max = 10 means we’ll consider p = 1..10
sel <- VARselect(data, lag.max = 10, type = "const")
p_aic <- sel$selection["AIC(n)"]
cat("VAR order by AIC:", p_aic, "\n")

# 2) Fit the VAR(p) with that lag
mod_var <- VAR(data, p = p_aic, type = "const")

# 3) Extract the residuals: this is a (T - p_aic) × 30 matrix
eps_mat <- residuals(mod_var)

acf(
  eps_mat[, 1],
  main = paste0("ACF of VAR(", p_aic, ") residuals,\nseries 1"),
  ylim = c(-1, 1)
)
Box.test(eps_mat[, 1], lag = 10, type = "Ljung-Box")


##############Fit the forest 

data <- eps_mat



lambda_grid <- seq(0.01, 0.2 ,  by = 0.01)
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
res0.4 <- main_oversteps_application( data = data , lambda_grid , grid = Grid_points_log,  start = NULL , 
                                      type = "SSR_row_log", k, p, num_class , cl )