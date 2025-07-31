install.packages("mev")
install.packages("tailDepFun")
install.packages("graphicalExtremes")
install.packages("igraph")
install.packages("mvtnorm")
library("mev")
library("tailDepFun")
library("graphicalExtremes")
library("igraph")
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


A <- rbind(c(1/3 , 1/3 , 1/3 , 0),
           c(1/2 , 0 , 0 , 1/2),
           c(0 , 3/4 , 1/4 , 0) ,
           c(0 , 1/2 , 0 , 1/2)
           )

# Sampling from the prior
# K: number of samples to draw from the prior


sampler <- function(K) {
  a11 <- runif(K , min = 0 , max = 1)
  a12 <- runif(K , min = 0 , max = 1)
  a21 <- runif(K , min = 0 , max = 1)
  a22 <-  runif(K , min = 0 , max = 1)
  alpha <- runif(K , min = 0 , max = 1)
  posterior <- matrix( c(alpha, a11 , a12 , a21 , a22) , byrow = TRUE, ncol = K)
  return(posterior)
}




# Marginal simulation from the statistical model
# theta: a matrix of parameters drawn from the prior
# m: number of conditionally independent replicates for each parameter vector


simulate <- function(Theta, m) {
  apply(Theta, 2, function(theta) {
    A <- matrix( theta[1:4] , byrow = TRUE , nrow = 2  )
    Z <- N_generate_Mix_log(m , A , theta[5])
    t(Z)
  }, simplify = FALSE)
}


m <- 30
Z_train <- simulate(theta_train, m)
Z_val   <- simulate(theta_val, m)

# Initialise the estimator
estimator <- juliaEval('

  d = 2    # dimension of each replicate
  w = 32   # number of neurons in each hidden layer

  # Final layer for 5 parameters in (0, 1)
  final_layer = Parallel(
    vcat,
    Dense(w, 1, σ),
    Dense(w, 1, σ),
    Dense(w, 1, σ),
    Dense(w, 1, σ),
    Dense(w, 1, σ)
  )

  psi = Chain(Dense(d, w, relu), Dense(w, w, relu), Dense(w, w, relu))
  phi = Chain(Dense(w, w, relu), Dense(w, w, relu), final_layer)
  deepset = DeepSet(psi, phi)
  estimator = PointEstimator(deepset)
')


K <- 5000
m <- 250
theta_train <- sampler(K)
theta_train[2,1:500] <- 0
theta_val   <- sampler(K/10)
theta_val[2,1:500] <- 0


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


theta_test <- matrix( c(0.25 , 0 , 1 , 1/2 , 1/2 ) , ncol = 1)
Z_test     <- simulate(theta_test, m)
assessment <- assess(estimator, theta_test, Z_test, estimator_names = "NBE",  parameter_names = c("alpha", "a11", "a12", "a21" , "a22"))
#>  Running NBE...
plotestimates(assessment)
