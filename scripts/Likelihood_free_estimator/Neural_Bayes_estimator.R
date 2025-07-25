library("NeuralEstimators")
library("JuliaConnectoR")
library("ggplot2")
library("ggpubr")
library("latex2exp")
juliaEval('using NeuralEstimators, Flux')


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

mu    <- theta_train[1, 1:4]
sigma <- theta_train[2, 1:4]
Z     <- Z_train[1:4]

df <- Map(function(z, m, s) {
  data.frame(Z = c(z), mu = m, sigma = s)
}, Z, mu, sigma)
df <- do.call(rbind, df)

df$theta <- paste0("$\\mu$ = ", round(df$mu, 3), ", $\\sigma$ = ", round(df$sigma, 3))
df$theta <- as.factor(df$theta)
levels(df$theta) <- sapply(levels(df$theta), TeX)

ggplot(df) +
  geom_histogram(aes(x = Z), bins = 30, fill = "blue", alpha = 0.5, color = "black") +
  facet_grid(~theta, labeller = label_parsed) +
  theme_bw()

# Initialise the estimator
estimator <- juliaEval('
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
  psi = Chain(Dense(n, w, relu), Dense(w, w, relu))
  phi = Chain(Dense(w, w, relu), final_layer)
  network = DeepSet(psi, phi)

  # Wrap the neural network as a PointEstimator
  estimator = PointEstimator(network)
')


# Train the estimator
estimator <- train(
  estimator,
  theta_train = theta_train,
  theta_val   = theta_val,
  Z_train = Z_train,
  Z_val   = Z_val
)


# Assess the estimator
theta_test <- sampler(10)
Z_test     <- simulate(theta_test, m)
assessment <- assess(estimator, theta_test, Z_test, estimator_names = "NBE")
#>  Running NBE...
plotestimates(assessment, parameter_labels = c("θ1" = expression(mu), "θ2" = expression(sigma)))
