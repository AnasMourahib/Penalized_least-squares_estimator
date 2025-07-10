#' K-Fold Cross-Validation Score for Penalized Estimation
#'
#' Computes the K=10 cross-validation score corresponding to a given penalization parameter \code{lambda}, based on Algorithm 2 from Mourahib, Kiriliouk, and Segers (2025). 
#' The function evaluates the model fit by repeatedly estimating the penalized least squares estimator on training folds and computing prediction error on the corresponding test folds.
#'
#' @param d Integer. Dimension of the problem (number of variables).
#' @param r Integer. Number of columns in the true coefficient matrix \code{A}.
#' @param A Numeric matrix. The true coefficient matrix.
#' @param grid Numeric matrix. Evaluation grid, where each row is a point \eqn{c_m}.
#' @param lambda Numeric. Penalization parameter for the estimation procedure.
#' @param num_col Integer (optional). Number of columns (i.e., number of extreme directions) to be used in the estimation. Defaults to \code{r}.
#' @param start Numeric vector. A \code{d × num_col} vectorized starting value for the optimization procedure. The rows of the coefficient matrix should be stacked sequentially.
#' @param type Character string. Specifies the fitted model type: either \code{"SSR_row_log"} for the logistic mixture model or \code{"SSR_row_HR"} for the Hüsler--Reiss mixture model.
#' @param p Numeric. Penalization exponent.
#' @param w A list containing two sub-lists: \code{w$train} and \code{w$test}. Each sub-list contains \code{num_class} elements (typically 10), each of which is a numeric vector of length \code{q} (number of evaluation points). For fold \code{i}, \code{w$train[[i]]} contains empirical stdf estimates for the q points based on the training set \eqn{S \setminus S_i}, while \code{w$test[[i]]} contains those based on the validation set \eqn{S_i}.
#' @param num_class Integer. Number of cross-validation folds. Default is 10.
#'
#' @return A numeric scalar: the mean cross-validation score computed across all folds, as described in Algorithm 2 from Mourahib, Anas, Anna Kiriliouk, and Johan Segers (2025), *A penalized least squares estimator for extreme-value mixture models*, arXiv preprint arXiv:2506.15272.
#'
#' @references
#' Mourahib, A., Kiriliouk, A., & Segers, J. (2025). A penalized least squares estimator for extreme-value mixture models. arXiv preprint arXiv:2506.15272.
#'
#' @export  


cross_validation<-function(d , r, A , grid, lambda, num_col = NULL, start , type = c("SSR_row_HR", "SSR_row_log"), p , w , num_class=10){
  print(lambda)
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

