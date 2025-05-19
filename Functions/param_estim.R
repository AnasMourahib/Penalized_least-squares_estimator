param_estim <- function(d , r , A , grid , lambda , num_col = NULL , start , type = c("SSR_row_HR", "SSR_row_log"), p ,  w ){
  q <- nrow(grid)
  if(is.null(num_col)){num_col <- r} #if not specified, use the correct number of columns
  l <- d*num_col
  v <- ifelse(type == "SSR_row_log", 1,  (d^(2) - d)/2   )
  interm <- function(theta, w){
    theta_A <- theta[1:l]
     if( any(  theta_A >1 | theta_A <0   ) ){
        return(10^16)
    }
    if(type=="SSR_row_log"){
      theta_alpha <- theta[l+v]
      if( theta_alpha>1 ) { return(10^16) }
      result <- .C( type ,  as.double(p),as.double(lambda), as.double(theta_A), as.integer(d), 
                    as.integer(num_col), as.integer(q),  as.double(rep(theta_alpha, num_col)), as.double(w), as.double(c(t(grid))), R = double(1))$R 
    }
    if(type=="SSR_row_HR"){
      theta_Gamma <-     theta[ (l+1) :  (l+ v)]
      if( length( which(theta_Gamma<0) ) > 0 ) { return(10^16) }
      result <- .C( type ,  as.double(p),as.double(lambda), as.double(theta_A), as.integer(d), 
                    as.integer(num_col), as.integer(q),  as.double(theta_Gamma), as.double(w), as.double(c(t(grid))), R = double(1))$R 
    }
    if (!is.finite(result)) {
      return(10^16)  # Return a large penalty value to handle non-finite values
    }
    return(result)
  }
  if(v == 1){
    start_dep <- runif(1 , min=0.1 , max = 0.9)
  } 
  else{
    start_dep <- runif( v  , min = 0.2 , max = 2)
  }
  start <- c(start , start_dep) 
  temp <- optim( start , interm , w = w , method="L-BFGS-B" , lower = rep(0, (l+v)  ) , control = list('maxit' = 1000))
  estim_pr <- temp$par
  theta_A_pr <- estim_pr[1:l]  #This is the estimation asscociated to the matrix
  A_pr <- matrix(   theta_A_pr , nrow = d , byrow = TRUE  )  
  theta_Z_pr <- estim_pr[(l+1) : (l+v)]
  ################
  
  interm_2 <- function(theta , w){
    if( any(  theta >1 | theta <0   ) ){
         return(10^16)
     }
    if (type == "SSR_row_log"){
      result <- .C(type , as.double(p),as.double(lambda), as.double(theta), as.integer(d), 
                   as.integer(num_col), as.integer(q), as.double(rep(dep, num_col)), as.double(w), as.double(c(t(grid))), R = double(1))$R  
    }
    if(type == "SSR_row_HR"){
      result <- .C(type , as.double(p),as.double(lambda), as.double(theta), as.integer(d), 
                   as.integer(num_col), as.integer(q), as.double(theta_Z_pr), as.double(w), as.double(c(t(grid))), R = double(1))$R  
    }
    if (!is.finite(result)) {
      return(10^16)  # Return a large penalty value to handle non-finite values
    }
    return(result)
  }
  
  start <- theta_A_pr
  temp <- optim( start , interm_2 , w = w , method="L-BFGS-B" , lower=rep(0 , l) , control = list('maxit' = 1000) )
  PLS_A <- temp$par
  matrix_A_PLS <- matrix( normalize_group( PLS_A , num_col )  , ncol = num_col , byrow = T)  
  if(type == "SSR_row_HR"){
    theta_Z_pr <- construct_symmetric_matrix(theta_Z_pr)
  }
  return(list("pls_matrix" =  matrix_A_PLS ,  "pls_dep" = theta_Z_pr ))
}
