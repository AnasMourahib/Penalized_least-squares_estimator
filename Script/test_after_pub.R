library("mev")
library("graphicalExtremes")
library("tailDepFun")
library("parallel")
library(mvtnorm)
dyn.load(normalizePath("main.dll"))
source("C:/Github/Penalized_least-squares_estimator/R_functions/helper_functions.R")
source("C:/Github/Penalized_least-squares_estimator/R_functions/estimation.R")

#####Evaluate EDS when sampling and fitting mixture logistic model
branch_path <- "C:/Results_Github/Penalized_least_squares_estimator/Extreme_directions_identificcation/DA_mix_log/Frechet3_Noise/"

n <- 100
N <- seq(500, 5000, by = 500)
lN <- length(N)

score_EDS <- numeric(lN)

for (i in 1:n) {
  file_path <- file.path(branch_path, paste0("result_", i, ".rds"))
  res <- readRDS(file_path)

  true_matrix <- res$True_matrixA
  true_extreme_directions <- extract_signatures(true_matrix)

  for (j in seq_along(N)) {  # seq_along(N)
    N_name <- paste0("N_", N[j])

    estim_matrix <- res$Estimation[[N_name]]$
      Estimation$kn_0.05$Estimation$pls_matrix

    estim_extreme_directions <- extract_signatures(estim_matrix)

    EDSij <- EDS(
      true_extreme_directions,
      estim_extreme_directions
    )
    cat("for see" , i , "this is the EDS" , EDSij ,"\n")
    score_EDS[j] <- score_EDS[j] + EDSij / n
  }
}

names(score_EDS) <- paste0("N_", N)
score_EDS





#####Estimation of a mixture logistic model
A <- rbind(
  c(1/3 , 1/3 , 0 , 1/3) ,
  c(1/2 , 1/2 , 0 , 0 ) ,
  c(1/2 , 0 , 0 , 1/2) ,
  c(1/3 , 0 , 1/3 , 1/3) ,
  c(1/2 , 0 , 1/2 , 0 )
)
Gamma <- rbind(
  c(0.0, 0.2, 0.4,  NA,  NA),
  c(0.2, 0.0, 0.6,  NA,  NA),
  c(0.4, 0.6, 0.0, 0.2, 0.6),
  c( NA,  NA, 0.2, 0.0, 0.4),
  c( NA,  NA, 0.6, 0.4, 0.0)
)
Gamma <- complete_Gamma(Gamma)
print(Gamma)
N <- 3000
Sigma <- sqrt(Gamma)/2
set.seed(7)
X <- N_generate_Mix_hr(N , A , Gamma)
X <- X[ , c(1 , 2 , 3)]
head(X)

#####Test main_oversteps_function





N <- nrow(X)
d <- ncol(X)

points <- c(0 , 1/6,  1/8 ,  1/4 ,  1/3 , 1/2 , 2/3 , 3/4 , 1)
grid <- selectGrid(points, d = d, nonzero = 2  )
q_HR <- nrow(grid)


R <- apply(X, 2L, rank)
kn <- 0.1
r <- 5L
p <- 0.4



num_class <- 5L

lambda_star_i <- 1 / sqrt(kn * N)

lambda_grid <- 1
print(lambda_grid)

stopifnot(
  ncol(grid) == d,
  nrow(R) == N,
  ncol(R) == d,
  k >= 1L,
  k < N
)

cl <- makeCluster(2L)

setup_estimation_workers(
  cl = cl,
  main_file = "R_functions/estimation.R",
  helper_file = "R_functions/helper_functions.R",
  dll_file = "main.dll"
)

clusterExport(cl, "stdfEmp", envir = .GlobalEnv)

fit_test <- main_oversteps(
  lambda_grid = lambda_grid,
  N = N,
  grid = grid,
  start = NULL,
  type = "SSR_row_HR",
  k = kn * N,
  p = p,
  num_class = num_class,
  cl = cl,
  d = d,
  r = r,
  sparse_proportion = NULL,
  X = X,
  q = q_HR,
  R = R,
  task = "ED_identification",
  seed = 1234L,
  max_num_col = 5L,
  maxit_cv = 50L,
  maxit_final = 100L,
  refined_grid = 0,
  verbose = TRUE
)

stopCluster(cl)

fit_test$num_col
fit_test$lambda_optim
fit_test$Estimation$pls_matrix
fit_test$Estimation$pls_dep
fit_test$cv_scores




#######This is only for analyzing the results, remove it before publication


d <- 3
points_HR <- c(0 , runif(1 , 0.2 , 0.4)   , runif(1 , 0.6 , 0.8) , runif(1 , 0.8 , 0.9) ,1)
Grid_points_HR <- selectGrid(cst = points_HR, d = d, nonzero = c(2,3))
head(Grid_points_HR)
q_HR <- nrow(Grid_points_HR)
print(q_HR)

q_old <- nrow(selectGrid(cst = points_HR, d = 3, nonzero = 2))
q_new <- nrow(selectGrid(cst = points_HR, d = 3, nonzero = c(2, 3)))

q_old
q_new
q_new / q_old


######Evaluation of each seed seperatly. Remove this before publication
true_matrix <- result_3$True_matrixA

estim_matrix <- result_3$Estimation[[2]]$Estimation$pls_matrix
estim_Gamma <-  result_3$Estimation[[2]]$Estimation$pls_dep

Gamma <- matrix(c(
  0,     0.08,  0.14,
  0.08,  0,     0.22,
  0.14,  0.22,  0
), nrow = 3, byrow = TRUE)


lambda <- result_1$Estimation[[2]]$lambda_min

lambda_grid <- result_1$Estimation[[2]]$lambda_grid



print(true_matrix)
print(estim_matrix)
print(Gamma)
print(estim_Gamma)
print(lambda_grid)
print(lambda)
