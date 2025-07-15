library(parallel)
library(graphicalExtremes)
library(mev) # generating of multivariate extreme-value data
library(gtools) # calculating permutations
library(tailDepFun) # defining a grid
library(ggplot2) # plotting
library(tidyverse)
dyn.load("main.dll")
source("R/Simulation_mixture_model.R")
source("R/Starting_points.R")
source("R/Help_functions.R")
source("R/param_estim.R")
source("R/cross_validation.R")
source("R/main.R")
source("R/main_applications.R")






#####Extracting data


X <- danube$data_clustered[ , c(7 , 18 , 24 , 27 , 30) ]
X <- as.matrix(X)
rownames(X) <- NULL
colnames(X) <- NULL
d <- ncol(X)
n <- nrow(X)



####Choice of tuning parameters

lambda_grid <- seq(0.01, 1 , length.out = 50)
p <- 0.4
k <- nrow(X)/10
num_class <- 5
points_log <- c(0,1/4 , 1/3 , 1/2 , 3/4 ,1)
Grid_points_log <- selectGrid(cst = points_log, d = d, nonzero  = c( 2,3) )

cl <- makeCluster(6)
# Export necessary functions and variables to cluster
clusterExport(cl, varlist = c(   "main_application", "shuffleCols" ,"param_estim_application", "construct_symmetric_matrix",  "normalize_group", "cross_validation_application", "p", "k","Grid_points_log" , "X"  , "num_class"
                                 , "lambda_grid"
))
clusterEvalQ(cl, { dyn.load("main.dll") })

set.seed(123)




#Detect the extreme directions
#res <- main_oversteps_application( data = X , lambda_grid , grid = Grid_points_log,  start = NULL ,
#                                      type = "SSR_row_log", k, p, num_class , cl )

###To directly get the results res <- readRDS(file = "Results/application/Danube_res0.4.rds")
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


#####Performance assessement
######## These are the estimated extremal correlations
mat_ext_corr_estim <- matrix(0 , d , d)
f<- function(vec) sum(vec^(1/coeff)  )^(coeff)
for(s in 1 : d){
  for(t in 1 : d){
    mat_ext_corr_estim[s,t] <- sum(  apply(mat[c(s,t) , ] ,  2 , f) )
  }
}

matrix_chi_estim <- 2 - mat_ext_corr_estim




########Empirical extremal correlation matrix
q <- 0.9
empirical_cdf <- function(x) {
  # rank each x[i] by how many values are ≤ it,
  # then divide by n to get the empirical CDF at x[i]
  rank(x, ties.method = "max") / length(x)
}




matrix_chi_emp <- matrix( 0 , nrow = d , ncol = d)



for( s in 1:d){
  vec_s <- X[,s]
  for(t in 1:d){
    vec_t <- X[,t]
    cdf_s <- empirical_cdf(vec_s)
    cdf_t <- empirical_cdf(vec_t)
    chi_st <- (length(which( (cdf_s >q)  &    (cdf_t >q) )    ) / n )   / (1 - q)
    if(chi_st >1) {chi_st <- 1}
    matrix_chi_emp[s,t] <- chi_st
  }
}

chi_estim <- matrix_chi_estim[upper.tri(matrix_chi_estim)]
chi_emp <- matrix_chi_emp[upper.tri(matrix_chi_emp)]

library(ggplot2)

df <- data.frame(
  fitted    = chi_estim,
  empirical = chi_emp
)

ggplot(df, aes(x = fitted, y = empirical)) +
  geom_point(
    color = "gray50",
    size  = 3,
    alpha = 0.8
  ) +
  geom_abline(
    intercept = 0, slope = 1,
    color     = "black",
    linetype  = "dashed",
    size      = 1
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "fitted", y = "empirical") +
  theme_minimal(base_size = 14) +      # base text size
  theme(
    panel.grid   = element_blank(),
    axis.line.x  = element_line(color = "black"),
    axis.line.y  = element_line(color = "black"),
    axis.title.x = element_text(size = 22),  # axis‑title font size
    axis.title.y = element_text(size = 22),
    axis.text.x  = element_text(size = 22),  # tick‑label font size
    axis.text.y  = element_text(size = 22)
  )


#########Estimate the probability of a failure set


res <- readRDS("C:/Users/mourahib/Desktop/github/PenalizedLeastSquaresEstimator/Results/application/Danube_res0.4.rds")
mat <- res$Estimation$matrix
alpha <- res$Estimation$dep

X_pareto <- apply(X, 2, function(col) {
  # Empirical ranks (largest = 1)
  r <- rank(-col, ties.method = "average")
  n <- length(col)
  # Transform to Pareto(1): X_i = n / rank_i
  n / r
})

Sum_per_observ <- apply(X_pareto ,  1 , sum  )
#q0.9 is the 0.9 quantile
threshold <-  seq(0.5 , 1 , by  = 0.05)
lthreshold <- length(threshold)
n <- 10^4
N <- 10^2
estim_pr <- rep(0 , lthreshold)
for (j in 1 : lthreshold){
  q <- quantile(Sum_per_observ , probs = threshold[j])
  for( i in 1 : N){
    set.seed(i)
    X <- N_generate_Mix_log(n , mat , alpha)
    X_pareto <- apply(X, 2, function(col) {
      # Empirical ranks (largest = 1)
      r <- rank(-col, ties.method = "average")
      n <- length(col)
      # Transform to Pareto(1): X_i = n / rank_i
      n / r
    })
    Sum_per_observ_i <- apply(X_pareto ,  1 , sum  )
    estim_pr[j] <- estim_pr[j] + ( (length(which(Sum_per_observ_i > q )) / n ) /(N)    )
  }
}






