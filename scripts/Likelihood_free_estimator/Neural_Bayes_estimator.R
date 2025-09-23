install.packages("mev")
install.packages("tailDepFun")
install.packages("graphicalExtremes")
install.packages("igraph")
install.packages("mvtnorm")
install.packages("JuliaConnectoR")
library("mev")
library("tailDepFun")
library("tibble")
library(mvtnorm)
library(parallel)
library(Matrix)
library(dplyr)
source("Functions/Simulation_mixture_model.R")

library("NeuralEstimators")
library("JuliaConnectoR")
library("ggplot2")
library("ggpubr")
library("latex2exp")
juliaEval('using NeuralEstimators, Flux')




# Sampling from the prior
# K: number of samples to draw from the prior


sampler <- function(K) {
  posterior <- runif(K , min = 0 , max = 1)
  # matrix( c(alpha , a11 ) , byrow = TRUE, ncol = K)
  return(posterior)
}




# Marginal simulation from the statistical model
# theta: a matrix of parameters drawn from the prior
# m: number of conditionally independent replicates for each parameter vector

alpha <- 0.25

simulate <- function(Theta, m) {
  lapply(Theta ,  function(theta) {
    A <- matrix( c(theta , 1 , 1/2 , 1/2) , byrow = TRUE , nrow = 2  )
    Z <- N_generate_Mix_log(m , A , alpha)
    t(Z)
  }
  )
}





# Initialise the estimator
estimator <- juliaEval('

  d = 2    # dimension of each replicate
  w = 32   # number of neurons in each hidden layer

  # Final layer for one parameters in [0, 1]
  final_layer =  Dense(w, 1, hardσ)

  psi = Chain(Dense(d, w, relu), Dense(w, w, relu), Dense(w, w, relu))
  phi = Chain(Dense(w, w, relu), Dense(w, w, relu), final_layer)
  deepset = DeepSet(psi, phi)
  estimator = PointEstimator(deepset)
')


K <- 5000
m <- 250
theta_train <- sampler(K)
theta_train[1:500] <- 0
theta_val   <- sampler(K/10)
theta_val[1:50] <- 0


Z_train <- simulate(theta_train, m)
Z_val   <- simulate(theta_val, m)

# Train the estimator
estimator <- train(
  estimator,
  theta_train = theta_train,
  theta_val   = theta_val,
  Z_train = Z_train,
  Z_val   = Z_val,
  epochs = 20
)

K_test <- 1000
theta_test <- sampler(K_test)
theta_test[1:50] <- 0
Z_test     <- simulate(theta_test, m)
assessment <- assess(estimator, theta_test, Z_test, estimator_names = "NBE",  parameter_names = c( "a"))
#>  Running NBE...
plotestimates(assessment)














##########Test with penalization term


juliaEval('include("C:/Users/20254817/Desktop/Githib/Penalized_least-squares_estimator/scripts/Likelihood_free_estimator/penalized_loss.jl")')




#### Test that the Julia function penalized_loss is working


juliaEval("isdefined(Main, :penalized_loss)")
juliaEval("isdefined(Main, :custom_loss)")
# should return TRUE

juliaEval("penalized_loss([0.1], [0.0] )")
# should return about 0.1001
juliaEval("custom_label(0.1 , 0.0)")
####

estimator <- train(
  estimator,
  theta_train = theta_train,
  theta_val   = theta_val,
  Z_train = Z_train,
  Z_val   = Z_val,
  epochs = 20,
   loss = "custom_loss"
)
