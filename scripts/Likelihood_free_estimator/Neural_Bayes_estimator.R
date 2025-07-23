# Install packages if not already installed
install.packages("JuliaCall")
install.packages("JuliaConnectoR")

# Load libraries
library(JuliaConnectoR)  # optional if you use JuliaCall primarily
library(JuliaCall)

# Setup Julia with the folder containing julia.exe
julia_setup(JULIA_HOME = "C:/Users/mourahib/AppData/Local/Programs/Julia-1.11.6/bin")

# Install packages in Julia if not already installed
julia_command('using Pkg; Pkg.add("NeuralEstimators")')
julia_command('using Pkg; Pkg.add("Flux")')

# Load the packages in the current Julia session
julia_library("NeuralEstimators")
julia_library("Flux")

# Alternatively, you can combine loading packages in one command
# juliaEval("using NeuralEstimators, Flux")



#> Starting Julia ...
#>



####Neural Bayes estimator for mean and standard deviation of a Gaussian random variable

# Sampling from the prior
# K: number of samples to draw from the prior
sampler <- function(K) {
  mu    <- rnorm(K)
  sigma <- rgamma(K, 1)
  Theta <- matrix(c(mu, sigma), byrow = TRUE, ncol = K)
  return(Theta)
}
K <- 10000
theta_train <- sampler(K)
theta_val   <- sampler(K/10)


# Marginal simulation from the statistical model
# theta: a matrix of parameters drawn from the prior
# m: number of conditionally independent replicates for each parameter vector
simulate <- function(Theta, m) {
  apply(Theta, 2, function(theta) t(rnorm(m, theta[1], theta[2])), simplify = FALSE)
}

m <- 30
Z_train <- simulate(theta_train, m)
Z_val   <- simulate(theta_val, m)


#Here we use a Neural network mapping from a simulation (m similations). Since our data eare replicated, we use Deepset structure to ensure
#excheangability. This deepset consists of two neural networks: an innter nn which has input of dimension d and outputs a dimension (d also probably). This neural network $
#gives summary statistics. For data withut spatial or temporal dependency, this neural network is a multilayer percpetron. The outer network maps
#the (m) learned summaries. This is always a multileyr perceptron.

# Initialise the estimator
estimator <- juliaEval('
  using NeuralEstimators: Parallel, DeepSet, PointEstimator
  using Flux: Dense, Chain, relu, softplus
  n = 1    # dimension of each data replicate (univariate)
  d = 2    # dimension of the parameter vector θ
  w = 128  # width of each hidden layer

  # Final layer has output dimension d and enforces parameter constraints
  final_layer = Parallel(
      vcat,
      Dense(w, 1, identity),     # mean parameter: no constraints
      Dense(w, 1, softplus)      # standard-deviation parameter: strictly positive
  )

  # Construct inner and outer networks and combine into DeepSet
  psi = Chain(Dense(n, w, relu), Dense(w, w, relu))  #This means that my inner neural network is a perceptron with 2 layers, each layer has 128 nodes and the activation funciton is the Relu
  phi = Chain(Dense(w, w, relu), final_layer)         #
  network = DeepSet(psi, phi)

  # Wrap the neural network as a PointEstimator
  estimator = PointEstimator(network)
')


juliaEval("using NeuralEstimators")
juliaEval("train")  # optional check — will throw an error if still unavailable


# Train the estimator
# assuming estimator, theta_train, Z_train, theta_val, Z_val are properly assigned in Julia:

trained_estimator <- julia_call(
  "train",
  estimator,
  theta_train,
  Z_train,
  theta_val = theta_val,
  Z_val = Z_val,
  epochs = as.integer(100),
  batchsize = as.integer(32),
  loss = "mse"  # or another valid loss function if required
)


