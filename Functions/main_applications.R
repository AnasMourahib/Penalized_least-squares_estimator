
param_estim_application <- function(d , grid , lambda , num_col , start , type = c("SSR_row_HR", "SSR_row_log"), p ,  w ){
  q <- nrow(grid)
  l <- d*num_col
  v <- numeric(l)
  #dim_dep is the degree of freedom of the variogram matrix 
  dim_dep <- 1
  #dim_dep <- 1
  print(dim_dep)
  interm <- function(theta, w){
    theta_A <- theta[1:l]
    theta_alpha <- theta[l+dim_dep]
    if( theta_alpha>1 ) { return(10^16) }
    result <- .C( type ,  as.double(p),as.double(lambda), as.double(theta_A), as.integer(d), 
                  as.integer(num_col), as.integer(q),  as.double(rep(theta_alpha, num_col)), as.double(w), as.double(c(t(grid))), R = double(1))$R 
    
    if (!is.finite(result)) {
      return(10^16)  # Return a large penalty value to handle non-finite values
    }
    return(result)
  }
  start_dep <- 0.5
  
  start <- c(start, start_dep)
  temp <- optim( start , interm , w = w , method="L-BFGS-B" , lower=rep(0, (l+dim_dep)  ) , control = list('maxit' = 1000))
  estim_fs <- temp$par
  #print(estim_fs)
  v_fs <- estim_fs[1:l]
  matrix <- matrix(   v_fs , nrow = d , byrow = TRUE  )  
  dep <- estim_fs[(l+1) : (l+dim_dep)]
  print(dep)
  ################
  
  interm_2 <- function(theta , w){
    
    result <- .C(type , as.double(p),as.double(lambda), as.double(theta), as.integer(d), 
                 as.integer(num_col), as.integer(q), as.double(rep(dep, num_col)), as.double(w), as.double(c(t(grid))), R = double(1))$R  
    
    if (!is.finite(result)) {
      return(10^16)  # Return a large penalty value to handle non-finite values
    }
    return(result)
  }
  
  start <- v_fs
  temp <- optim( start , interm_2 , w = w , method="L-BFGS-B" , lower=rep(0 , l) , control = list('maxit' = 1000) )
  estim <- temp$par
  matrix_estim <- matrix( normalize_group( estim , num_col )  , ncol = num_col , byrow = T) 
  
  
  
  
  return(list("matrix" = matrix_estim ,  "dep" = dep ))
}


cross_validation_application<-function(d, grid, lambda , num_col ,  start , type =  "SSR_row_log", p , w , num_class=5){
  #start <- c(start, 0.5)    
  w_train <- w$train 
  w_test <- w$test
  q <- nrow(grid)
  l <- d * num_col 
  CV <- vector(length = num_class)
  for (class_k in 1:num_class){
    optimizer_minus_class_k <- param_estim_application(d , grid , lambda , num_col , start , type = type , p ,  w = w_train[[class_k]])
    print(optimizer_minus_class_k$matrix)
    v_A <- as.vector( t(optimizer_minus_class_k$matrix) )
    v_alpha <- optimizer_minus_class_k$dep 
    
    CV[class_k] <- .C( type , as.double(p) , as.double(0) , as.double(v_A) , as.integer(d), 
                       as.integer(num_col) , as.integer(q) , as.double(rep(v_alpha, num_col)) , as.double(w_test[[class_k]]) , as.double(c(t(grid))) , R = double(1))$R 
    
  }
  return(mean(CV))
}





main_application <- function( data , lambda_grid , grid  , num_col , start = NULL,
                              type = "SSR_row_log" , k, p, num_class = num_class, cl) {
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
                                        type = "SSR_row_log", k, p, num_class , cl ){
  num_col <- 1
  
  Estimation <- main_application( data , lambda_grid , grid  , num_col , start = start,
                                  type = "SSR_row_log" , k, p, num_class = num_class, cl)
  test <- TRUE
  while(test == TRUE){
    #print(num_col)
    old_Estimation <- Estimation
    old_matrix <- Estimation$Estimation$matrix
    print(old_Estimation)
    sig_old <- apply(old_matrix, 2,  function(vec)   which(vec>0), simplify = FALSE )
    print(paste("The old estimation with", num_col, " number of columns is "))
    print(sig_old)
    #print(old_Estimation)
    num_col <- num_col + 1 
    Estimation <- main_application( data , lambda_grid , grid  , num_col , start = start,
                                    type = "SSR_row_log" , k, p, num_class = num_class, cl)
    estim_matrix <- Estimation$Estimation$matrix
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




#############################