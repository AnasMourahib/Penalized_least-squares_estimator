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
                                        type = "SSR_row_log", k, p, num_class , cl ){
  num_col <- 1
  
  Estimation <- main_application( data , lambda_grid , grid  , num_col , start = start,
                                  type = "SSR_row_log" , k, p, num_class = num_class, cl)
  test <- TRUE
  while(test == TRUE){
    #print(num_col)
    old_Estimation <- Estimation
    old_matrix <- Estimation$Estimation$matrix
    sig_old <- apply(old_matrix, 2,  function(vec)   which(vec>0), simplify = FALSE )
    print(paste("The old estimation with", num_col, " number of columns is "))
    print(sig_old)
    print(old_Estimation)
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