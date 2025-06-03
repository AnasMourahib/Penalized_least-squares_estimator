param_estim_application <- function(d, grid, lambda, num_col , start, type = c("SSR_row_HR", "SSR_row_log"), p, w) {
  q <- nrow(grid)
  l <- d * num_col
  v <- ifelse(type == "SSR_row_log", 1, (d * (d - 1)) / 2)
  
  grid_flat <- as.double(c(t(grid)))  # precompute
  lower_bounds <- rep(0, l + v)
  
  # interm function
  interm <- function(theta) {
    theta_A <- theta[1:l]
    if (any(theta_A > 1 | theta_A < 0)) return(1e16)
    
    if (type == "SSR_row_log") {
      theta_alpha <- theta[l + 1]
      if (theta_alpha > 1) return(1e16)
      result <- .C(type, as.double(p), as.double(lambda), as.double(theta_A), as.integer(d),
                   as.integer(num_col), as.integer(q), as.double(rep(theta_alpha, num_col)), as.double(w),
                   grid_flat, R = double(1))$R
    } else {
      theta_Gamma <- theta[(l + 1):(l + v)]
      if (any(theta_Gamma < 0)) return(1e16)
      result <- .C(type, as.double(p), as.double(lambda), as.double(theta_A), as.integer(d),
                   as.integer(num_col), as.integer(q), as.double(theta_Gamma), as.double(w),
                   grid_flat, R = double(1))$R
    }
    
    if (!is.finite(result)) return(1e16)
    result
  }
  
  # Initialization
  start_dep <- if (v == 1) runif(1, 0.1, 0.9) else rep(1 , v)
  start_total <- c(start, start_dep)
  
  temp <- optim(start_total, interm, method = "L-BFGS-B", lower = lower_bounds, control = list(maxit = 1000))
  estim_pr <- temp$par
  theta_A_pr <- estim_pr[1:l]
  A_pr <- matrix(theta_A_pr, nrow = d, byrow = TRUE)
  theta_Z_pr <- estim_pr[(l + 1):(l + v)]
  # Reuse variables
  A_vec <- as.vector(t(matrix(normalize_group(theta_A_pr, num_col), ncol = num_col, byrow = TRUE)))
  
  # Partial Estimation of A
  interm_2 <- function(theta) {
    if (any(theta > 1 | theta < 0)) return(1e16)
    dep <- if (type == "SSR_row_log") rep(1, num_col) else theta_Z_pr
    
    result <- .C(type, as.double(p), as.double(lambda), as.double(theta), as.integer(d),
                 as.integer(num_col), as.integer(q), as.double(dep), as.double(w),
                 grid_flat, R = double(1))$R
    
    if (!is.finite(result)) return(1e16)
    result
  }
  
  temp <- optim(theta_A_pr, interm_2, method = "L-BFGS-B", lower = rep(0, l), control = list(maxit = 1000))
  PLS_A <- temp$par
  matrix_A_PLS <- matrix(normalize_group(PLS_A, num_col), ncol = num_col, byrow = TRUE)
  return(list(pls_matrix = matrix_A_PLS, pls_dep = theta_Z_pr))
}



cross_validation_application<-function(d, grid, lambda , num_col ,  start , type =  type, p , w , num_class=5){
  #start <- c(start, 0.5)    
  w_train <- w$train 
  w_test <- w$test
  q <- nrow(grid)
  l <- d * num_col 
  CV <- vector(length = num_class)
  for (class_k in 1:num_class){
    optimizer_minus_class_k <- param_estim_application(d , grid , lambda , num_col , start , type = type , p ,  w = w_train[[class_k]])
    print(optimizer_minus_class_k$pls_matrix)
    v_A <- as.vector( t(optimizer_minus_class_k$pls_matrix) )
    v_dep <- optimizer_minus_class_k$pls_dep
    if(type == "SSR_row_log"){
      v_dep <- as.double(rep(v_dep  , num_col))
    }
    else{
      v_dep <- as.double(v_dep)
    }
    
    CV[class_k] <- .C( type , as.double(p) , as.double(0) , as.double(v_A) , as.integer(d), 
                       as.integer(num_col) , as.integer(q) , v_dep , as.double(w_test[[class_k]]) , as.double(c(t(grid))) , R = double(1))$R 
    
  }
  return(mean(CV))
}


main_application <- function( data , lambda_grid , grid  , num_col , start = NULL,
                              type = type , k, p, num_class = num_class, cl) {
  d <- ncol(data)  
  q <- nrow(grid)
  l <- d * num_col
  #checkN <- N %% num_class
  #if (checkN > 0) {
  #  N <- N - checkN
  #  cat("Warning: sample size N has been reduced to", N, "\n")
  #}
  
  # Data generation and preparation
  
  R <- apply(data, 2, rank)
  #Calculate the empirical stdf estimates on the grid points for the cross validation part
  w <- W_calculus(k = k, num_class = num_class, X = data, grid = grid, q = q)
  #Calculate the empitical stdf estimates on the grid points for the estimation part
  w_total <- sapply(1:q, function(m) stdfEmp(R, k, grid[m, ]))
  
  if (is.null(start)) {
    start <- starting_point(data, num_col)
    start <- c(t(start))
  }
  
  # Cross-validation wrapper
  wrapper_cross_validation <- function(lambda) {
    cross_validation_application(d, grid, lambda , num_col ,  start , type = type , p , w , num_class=num_class)
  }     
  
  # Parallelize cross-validation
  scores <- unlist(parLapply(cl, lambda_grid, wrapper_cross_validation))
  #print(scores)
  # Determine optimal lambda
  index <- which(scores == min(scores))
  lambda_optim <- lambda_grid[index]
  paste0("This is the vallue of the optiam lambda" , lambda_optim)
  
  # Wrapper for param_estim
  wrapper_param_estim <- function(lambda) {
    param_estim_application(d , grid , lambda , num_col , start , type = type , p ,  w_total )
  }  
  
  # Execute param_estim in parallel
  Estimation <- parLapply(cl, as.list(lambda_optim), wrapper_param_estim)[[1]]
  # Compile results
  result <- list(lambda_optim = lambda_optim, Estimation = Estimation )
  
  return(result)
}



main_oversteps_application <- function( data , lambda_grid, grid,  start , 
                                        type , k, p, num_class , cl ){
  num_col <- 1
  
  Estimation <- main_application( data , lambda_grid , grid  , num_col , start = start,
                                  type = type , k, p, num_class = num_class, cl)
  test <- TRUE
  while(test == TRUE){
    #print(num_col)
    old_Estimation <- Estimation
    old_matrix <- Estimation$Estimation$pls_matrix
    sig_old <- apply(old_matrix, 2,  function(vec)   which(vec>0), simplify = FALSE )
    print(paste("The old estimation with", num_col, " number of columns is "))
    print(sig_old)
    print(old_Estimation)
    num_col <- num_col + 1 
    Estimation <- main_application( data , lambda_grid , grid  , num_col , start = start,
                                    type = type , k, p, num_class = num_class, cl)
    estim_matrix <- Estimation$Estimation$pls_matrix
    sig_new  <- apply( estim_matrix, 2,  function(vec)   which(vec>0), simplify = FALSE )
    # print(sig_new)
    #print(paste("The new estimation with", num_col, " number of columns is ", sig_new))
    #if( has_duplicates(sig_new)  ) 
    #{ print(test)
    #test <- FALSE} 
    test <- all( apply( estim_matrix , 2 , function(vec) all(vec==0) ) == FALSE ) 
    print(test)
    
  }
  
  return(old_Estimation)
  
}


