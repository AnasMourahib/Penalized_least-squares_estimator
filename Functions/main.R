main <- function(lambda_grid, N, seed, A, grid, alpha = NULL, Sigma = NULL ,  num_col = NULL, start = NULL,
                 type = c("SSR_row_HR","SSR_row_log"), k, p, num_class = 10, cl) {
  r <- ncol(A)
  d <- nrow(A)   
  q <- nrow(grid)
  
  if (is.null(num_col)) { 
    num_col <- r  # Default to the correct number of columns
  }
  
  l <- d * num_col
  checkN <- N %% num_class
  if (checkN > 0) {
    N <- N - checkN
    cat("Warning: sample size N has been reduced to", N, "\n")
  }
  
  # Data generation and preparation
  set.seed(seed)
  if(type=="SSR_row_HR"){
    X <- N_generate_Mix_hr(N, A, Sigma)
  }
  if(type=="SSR_row_log"){
    X <- N_generate_Mix_log(N, A, alpha)
    print(dim(X))
  }
  
  ##Introduce the noise
  
  Noise <- matrix( rnorm(N * r , mean = 0 , sd = 3)  ,  nrow = N , ncol = r )  
  X <- X + Noise
  
  ####
  
  R <- apply(X, 2, rank)
  #Calculate the empirical stdf estimates on the grid points for the cross validation part
  w <- W_calculus(k = k, num_class = num_class, X = X, grid = grid, q = q)
  #Calculate the empitical stdf estimates on the grid points for the estimation part
  w_total <- sapply(1:q, function(m) stdfEmp(R, k, grid[m, ]))
  
  if (is.null(start)) {
    start <- starting_point(X, num_col)
    start <- c(t(start))
  }
  
  # Cross-validation wrapper
  wrapper_cross_validation <- function(lambda) {
    cross_validation(d = d, r = r, A = A , grid = grid, lambda = lambda, num_col = num_col, 
                       start = start, type = type, p = p, w = w, num_class = num_class)
  }     
  
  # Parallelize cross-validation
  scores <- unlist(parLapply(cl, lambda_grid, wrapper_cross_validation))
  
  
  
  
  nan <-  which(is.nan(scores)  ) 
  if(length(nan)>0){
    scores <- scores[ -nan   ]
  }
  # Determine optimal lambda
  index <- which(scores == min(scores))
  lambda_optim <- lambda_grid[index]
  
  
  
  # Wrapper for param_estim
  wrapper_param_estim <- function(lambda) {
    param_estim(d = d, r = r, A = A, grid = grid, lambda = lambda, num_col = num_col, 
                  start = start, type = type, p = p, w = w_total)
  }  
  
  # Execute param_estim in parallel
  Estimation <- parLapply(cl, as.list(lambda_optim), wrapper_param_estim)[[1]]
  # Compile results
  result <- list(lambda_optim = lambda_optim, Estimation = Estimation )
  
  return(result)
}

