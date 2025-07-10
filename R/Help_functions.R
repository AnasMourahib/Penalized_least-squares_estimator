
# Puts the columns of start in the "right" order (i.e. the order corresponding to A)
# works also if ncol(start) > ncol(A)
shuffleCols <- function(start, A){
  r <- ncol(A)
  k <- ncol(start)
  perms <- permutations(n = k, r = r)
  temp <- apply(perms, 1, function(j) sum( (start[,j] - A)^(2)  ))
  indx <- which(temp == min(temp))
  firstcols <- perms[indx,]
  lastcols <- setdiff(c(1:k),firstcols)
  return(start[,c(firstcols,lastcols)])
}


construct_symmetric_matrix <- function(vec) {
  # Determine d from the length of vec
  d <- (1 + sqrt(1 + 8 * length(vec))) / 2
  if (d != floor(d)) stop("Vector length is incorrect for a symmetric matrix")
  d <- as.integer(d)

  # Initialize d x d matrix with 1s on the diagonal
  mat <- diag(1, d, d)

  # Fill the upper triangle row by row
  index <- 1
  for (i in 1:(d-1)) {
    for (j in (i+1):d) {
      mat[i, j] <- vec[index]
      index <- index + 1
    }
  }

  # Make the matrix symmetric
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]

  return(mat)
}


normalize_group <- function(v, group_size) {
  v_normalized <- sapply(seq(1, length(v), by = group_size), function(i) {
    group <- v[i:(i + group_size - 1)]
    group_sum <- sum(group)

    if (group_sum == 0) {
      return(group)  # Return the original group if all values are zero
    } else {
      return(group / group_sum)  # Otherwise, normalize
    }
  })

  return(as.vector(v_normalized))
}




W_calculus<-function(k, num_class, X, grid, q){
  W_train<-list()
  W_test<-list()
  N<-nrow(X)
  for(class_k in 1: num_class) {
    low <- ((class_k-1)*(N/num_class)+1)
    up <- (class_k*(N/num_class))
    X_train=X[-c(low:up),]
    X_test=X[c(low:up),]
    R_train<-apply(X_train,2,rank)
    R_test<-apply(X_test,2,rank)
    W_train[[class_k]] <- sapply(c(1:q), function(m) stdfEmp(R_train, k = (1-1/num_class)*k, grid[m,]))
    W_test[[class_k]]  <- sapply(c(1:q), function(m) stdfEmp(R_test, k = (1/num_class)*k, grid[m,]))
  }
  return(list("train"=W_train, "test"=W_test))
}







































extr_coeff_mix_log <- function( A , alpha){
  Apower <- A^(alpha)
  sum_cols <-  apply(A , 2, sum)
  res <- sum(  sum_cols^(1/alpha) )
  return(res)
}

ext_coeff_HR <- function(Gamma){
  cdf1 <- pnorm(  Gamma/2 , sd= sqrt(Gamma) )
  cdf2 <- pnorm(  Gamma/2 , sd= sqrt(Gamma) )
  result <- cdf1 + cdf2
  return(result)
}


EDS <- function(B, A){
  sig_A <- apply(A, 2,  function(vec)   which(vec>0) )
  sig_B <- apply(B, 2, function(vec)   which(vec>0) , simplify = FALSE )
  return( 1 -   (   length(intersect(sig_A , sig_B)) / length(union(sig_A , sig_B))   )      )
}


SMSE_log <- function (list, A , dep){
  d <- nrow(A)
  r <- ncol(A)
  N <- length(list)
  list_matrices <- list()
  list_dep <- rep(0 , N)
  for(i in 1: N){
    list_matrices[[i]] <-  list[[i]]$matrix
    r_diff_estim <-  r - ncol(list_matrices[[i]])
    if(r_diff_estim > 0 ){
      for(e in 1 : r_diff_estim){
        list_matrices[[i]] <- cbind(list_matrices[[i]] , rep(0 , d)  )
      }
    }
    list_dep[i]  <- list[[i]]$dep
  }
  rmse <- 0
  array_stack <- simplify2array(list_matrices)
  max_matrix <- apply(array_stack, c(1, 2), max)
  max_matrix <- apply(max_matrix, c(1,2), f<-function(x) {
    if(x==0)
    {x <- 1}
    return(x) } ) #Correct when sd=0. This case happens only when a true zero is estimated to 0 in all seeds
    true_ext_coeff <-  ext_coeff_mix_log(A , dep)
    estim_ext_coeff <- sapply(list_dep ,  function(v) ext_coeff_mix_log(A  , v) )
    max_ext_coeff <- max(estim_ext_coeff)^2
    ext_coeff_estim_diff <- sum( ( estim_ext_coeff - true_ext_coeff)^(2)  )



  for( i in 1:N){
    diff <-   as.vector( (list_matrices[[i]] - A)^2 )  / (N * max_matrix)
    rmse <- rmse + ( sum(diff)  )
  }
  rmse <- rmse + (r *  ext_coeff_estim_diff)/(N  * max_ext_coeff)
  return(rmse)
}


SMSE_log_2 <- function (list, A , dep){
  d <- nrow(A)
  r <- ncol(A)
  N <- length(list)
  list_matrices <- list()
  list_dep <- rep(0 , N)
  for(i in 1: N){
    list_matrices[[i]] <-  list[[i]]$Estimation$matrix
    r_diff_estim <-  r - ncol(list_matrices[[i]])
    if(r_diff_estim > 0 ){
      for(e in 1 : r_diff_estim){
        list_matrices[[i]] <- cbind(list_matrices[[i]] , rep(0 , d)  )
      }
    }
    list_dep[i]  <- list[[i]]$Estimation$dep
  }
  rmse <- 0
  array_stack <- simplify2array(list_matrices)
  max_matrix <- apply(array_stack, c(1, 2), max)
  max_matrix <- apply(max_matrix, c(1,2), f<-function(x) {
    if(x==0)
    {x <- 1}
    return(x) } ) #Correct when sd=0. This case happens only when a true zero is estimated to 0 in all seeds
  max_dep_coeff <- max(list_dep)^2
  dep_diff <- sum( ( list_dep - dep)^(2)  )



  for( i in 1:N){
    diff <-   as.vector( (list_matrices[[i]] - A)^2 )  / (N * max_matrix)
    rmse <- rmse + ( sum(diff)  )
  }
  rmse <- rmse + (r * dep_diff)/(N  * max_dep_coeff)
  return(rmse)
}



ext_coeff_HR <- function(Gamma){
  cdf1 <- pnorm(  Gamma/2 , sd= sqrt(Gamma) )
  cdf2 <- pnorm(  Gamma/2 , sd= sqrt(Gamma) )
  result <- cdf1 + cdf2
  return(result)
}

matrix_ext <- function(mat){
  upper_values <- mat[upper.tri(mat, diag = FALSE)]
  values <- ext_coeff_HR(upper_values)
  return(values)
}

matrix_vector <- function(mat, ind = NULL){
  upper_values <- mat[upper.tri(mat, diag = FALSE)]
  if(!is.null(ind)){
    upper_values <- upper_values[-ind]
  }
  return(upper_values)
}

SMSE_HR <- function (list_matrices , A , Gamma  ){
  rmse <-  0
  d <- nrow(A)
  r <- ncol(A)
  N <- length(list_matrices)
  list_pls_matrices <- list()
  list_pls_dep_matrices <- list()
  for(i in 1: N){
    list_pls_matrices[[i]] <-  list_matrices[[i]]$Estimation$pls_matrix
    list_pls_dep_matrices[[i]] <- list_matrices[[i]]$Estimation$pls_dep
  }
  array_stack <- simplify2array(list_pls_matrices)
  max_matrix <- apply(array_stack, c(1, 2), max)
  max_matrix <- apply(max_matrix, c(1,2), f<-function(x) {
    if(x==0)
    {x <- 1}
    return(x) } ) #Correct when sd=0. This case happens only when a true zero is estimated to 0 in all seeds
    list_values <- lapply(list_pls_dep_matrices, matrix_ext)
    max_vector <- do.call(pmax, list_values)
    true_value <- ext_coeff_HR(Gamma)
    for( i in 1:N){
      diff1 <- sum( as.vector( (list_pls_matrices[[i]] - A)^2 )  / (N * max_matrix) )
      diff2 <- sum((list_values[[i]] - true_value)^2  / (N * max_vector))
      rmse <- rmse + (diff1 + diff2)
    }
    return(rmse)
}

SMSE_HR2 <- function (list_matrices , A , ind = NULL  ){
  rmse <-  0
  d <- nrow(A)
  r <- ncol(A)
  N <- length(list_matrices)
  list_pls_matrices <- list()
  list_pls_dep_matrices <- list()
  for(i in 1: N){
    list_pls_matrices[[i]] <-  list_matrices[[i]]$Estimation$pls_matrix
    list_pls_dep_matrices[[i]] <- list_matrices[[i]]$Estimation$pls_dep
  }
  array_stack <- simplify2array(list_pls_matrices)
  max_matrix <- apply(array_stack, c(1, 2), max)
  max_matrix <- apply(max_matrix, c(1,2), f<-function(x) {
    if(x==0)
    {x <- 1}
    return(x) } ) #Correct when sd=0. This case happens only when a true zero is estimated to 0 in all seeds
  interm_matrix_vector <- function(matrix){
    res <- matrix_vector(matrix , ind)
    return(res)
  }
  list_values <- lapply(list_pls_dep_matrices, interm_matrix_vector  )
  max_vector <- (do.call(pmax, list_values) )^2
  true_value <- 1
  for( i in 1:N){
    diff1 <- sum( as.vector( (list_pls_matrices[[i]] - A)^2 )  / (N * max_matrix) )
    diff2 <- sum((list_values[[i]] - true_value)^2  / (N * max_vector))
    rmse <- rmse + (diff1 + diff2)
  }
  return(rmse)
}

