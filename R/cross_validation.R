cross_validation<-function(d , r, A , grid, lambda, num_col = NULL, start , type = c("SSR_row_HR", "SSR_row_log"), p , w , num_class=10){
  print(lambda)
  #start <- c(start, 0.5)    
  w_train <- w$train 
  w_test <- w$test
  q <- nrow(grid)
  if(is.null(num_col)){num_col <- r} #if not specified, use the correct number of columns
  l <- d * num_col 
  CV <- vector(length = num_class)
  for (class_k in 1:num_class){
    optimizer_minus_class_k <- param_estim(d = d , r = r , A = A , grid = grid , lambda = lambda , num_col = num_col , start = start  , type = type , p = p  ,  w = w_train[[class_k]] )
    
    v_A <- as.vector( t(optimizer_minus_class_k$pls_matrix) )
    v_alpha <- optimizer_minus_class_k$pls_dep 
    if(type == "SSR_row_log"){
      CV[class_k] <- .C( type , as.double(p) , as.double(0) , as.double(v_A) , as.integer(d), 
                         as.integer(num_col) , as.integer(q) , as.double(rep(v_alpha, num_col)) , as.double(w_test[[class_k]]) , as.double(c(t(grid))) , R = double(1))$R 
    }
    if(type == "SSR_row_HR"){
      CV[class_k] <- .C( type , as.double(p) , as.double(0) , as.double(v_A) , as.integer(d), 
                         as.integer(num_col) , as.integer(q) , as.double(v_alpha) , as.double(w_test[[class_k]]) , as.double(c(t(grid))) , R = double(1))$R 
    }
  }
  return(mean(CV))
}

