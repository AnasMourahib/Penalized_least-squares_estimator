install.packages("mev")
install.packages("reticulate")
library(reticulate)
library(mev)
source_python("scripts/DAMEX_function.py")
source("Functions/Simulation_mixture_model.R")
source("Functions/Help_functions.R")



seed = 1:100
lseed <- length(seed)
A <- rbind(c(1/3,0,1/3,1/3), c(0,1/2,1/2,0), c(3/4,0,0,1/4), c(1/2,1/2,0,0))
r = ncol(A)
######## dependence coefficient for the mixture logistic model
alpha = 0.25
####### dependence coefficient for the mixture HR model
d = nrow(A)
lambda_sq <- 1/4
Sigma <- matrix(lambda_sq , nrow = d , ncol = d )
diag(Sigma) <- 0
######
mu = seq(0.0001 , 0.01 , by = 0.001)
lext_dir = length(mu)
size = seq(500  , 5000 , by = 500)
lsize = length(size)
ED_S <- rep(0 , lsize)

for (i in 1:lsize) {
  res <- 0
  N = size[i]
  for (j in seed){
    set.seed(j)
    ###Simulation from the mixture logistic model
    #X <- N_generate_Mix_log(N , A , alpha)
    ###Simulation from the mixture HR model
    X <- N_generate_Mix_hr(N, A, Sigma)
    Noise <- matrix( rnorm(N * r , mean = 0 , sd = 3)  ,  nrow = N , ncol = r )
    X <- X + Noise
    ext_dir = DAMEX_test(data = X  , mu = mu )
    EDS <- rep(0 , lext_dir)
    for (k  in 1:lext_dir) {
      EDS[k] = EDS_python(ext_dir[[k]] , A)
    }
    min_EDS = min(EDS)
    res <- res + (min_EDS/lseed)
  }

  ED_S[i] = res


}



saveRDS(ED_S , file = "C:/Users/20254817/Desktop/Githib/Penalized_least-squares_estimator/Results/EDS_DAMEX_mixlog.RDS")
saveRDS(ED_S , file = "C:/Users/20254817/Desktop/Githib/Penalized_least-squares_estimator/Results/EDS_DAMEX_mixHR.RDS")




#################### plot EDS for mix log and mix HR

