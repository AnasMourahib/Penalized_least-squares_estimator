#' Reorder Columns to Match a Target Matrix
#'
#' Reorders the columns of the input matrix \code{start} to best match the columns of the matrix \code{A}, based on minimizing the squared Euclidean distance without normalization, as in Equation (4.5) from Mourahib, Kiriliouk, and Segers (2025).
#'
#' @param start A numeric matrix with \code{d} rows. Typically an estimated coefficient matrix.
#' @param A A numeric matrix representing the true coefficient matrix.
#'
#' @return A reordered version of \code{start}, where the first \code{ncol(A)} columns are matched to \code{A} to minimize squared differences, and remaining columns are placed at the end.
#'
#' @references
#' Mourahib, A., Kiriliouk, A., & Segers, J. (2025). A penalized least squares estimator for extreme-value mixture models. \emph{arXiv preprint} arXiv:2506.15272.
#'
#' @export


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

#' Construct a Symmetric Matrix from a Vector, typicaly, we use this to contruct the variogram matrix based on the upper triangular entries.
#'
#' Converts a vector containing the upper triangular entries (excluding the diagonal) of a symmetric matrix into the full symmetric matrix with ones on the diagonal.
#'
#' @param vec A numeric vector of length \eqn{d(d - 1)/2}, representing the upper triangle (excluding the diagonal) of a \eqn{d \times d} symmetric matrix.
#'
#' @return A \eqn{d \times d} symmetric matrix with ones on the diagonal and the upper/lower triangular entries filled from \code{vec}.
#'
#' @export

construct_symmetric_matrix <- function(vec) {
  # Determine d from the length of vec
  d <- (1 + sqrt(1 + 8 * length(vec))) / 2
  if (d != floor(d)) stop("Vector length is incorrect for a symmetric matrix")
  d <- as.integer(d)

  # Initialize d x d matrix with 0s on the diagonal
  mat <- diag(0, d, d)

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


#' Normalize Elements of a Vector by Groups
#'
#' Splits a vector into groups of equal size and normalizes each group so that its elements sum to 1. If a group's sum is zero, it is left unchanged.
#'
#' @param v A numeric vector to be normalized.
#' @param group_size An integer specifying the size of each group.
#'
#' @return A numeric vector of the same length as \code{v}, where each group of \code{group_size} elements is normalized to sum to 1 (unless the group sum is zero).
#'
#' @export

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


#' Compute Empirical Stable Tail Dependence Functions for Cross-Validation
#'
#' Splits the dataset into \code{num_class} folds and computes the empirical stable tail dependence function
#' (STDF) estimates using the \code{stdfEmp} function for each training and testing set. Typically used
#' in a K-fold cross-validation setting to prepare input for penalized model fitting.
#'
#' @param k The order statistic threshold, typically set to \code{floor(n^\alpha)} where \code{n} is sample size.
#' @param num_class Number of folds (classes) to split the data into for cross-validation (e.g., 10).
#' @param X A numeric matrix of dimension \code{n x d}, containing the observed multivariate data.
#' @param grid A matrix of evaluation points \eqn{c_1, \ldots, c_q}, where each row corresponds to a point in the unit simplex.
#' @param q The number of evaluation points (i.e., the number of rows in \code{grid}).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{train}}{A list of length \code{num_class}, each containing a vector of length \code{q} with STDF estimates using training data.}
#'   \item{\code{test}}{A list of length \code{num_class}, each containing a vector of length \code{q} with STDF estimates using test data.}
#' }
#'
#' @seealso \code{\link{stdfEmp}} for computing empirical STDF values.
#'
#' @export



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



#' Compute the extreme directions score
#'
#' Calculates a Jaccard distance between matrices \code{B} and \code{A} based on the proportion of the non-similar signatures
#'
#' @param B A numeric matrix.
#' @param A A numeric matrix with the same number of columns as \code{B}.
#'
#' @return A numeric value between 0 and 1 representing the proportion of missidentified signatures (or equivalently extreme directions); See Equation (4.3) from Mourahib, Kiriliouk, and Segers (2025).
#'
#' @references
#' Mourahib, A., Kiriliouk, A., & Segers, J. (2025). A penalized least squares estimator for extreme-value mixture models. \emph{arXiv preprint} arXiv:2506.15272.
#'
#' @export
#'
EDS <- function(B, A){
  sig_A <- apply(A, 2,  function(vec)   which(vec>0) )
  sig_B <- apply(B, 2, function(vec)   which(vec>0) , simplify = FALSE )
  return( 1 -   (   length(intersect(sig_A , sig_B)) / length(union(sig_A , sig_B))   )      )
}




#' Standardized mean squarred error in Equation (4.4) from Mourahib, Kiriliouk, and Segers (2025) for the mixture logistic model
#'
#' Computes a standardized mean squared error (SMSE) between a list of estimated coefficient matrices and a true coefficient matrix,
#' adjusted for differences in the number of columns. It also incorporates the squared error of an associated dependence parameter.
#'
#' @param list A list of estimation results, where each element is expected to contain a component \code{Estimation$matrix} (estimated matrix)
#' and \code{Estimation$dep} (estimated dependence parameter).
#' @param A The true coefficient matrix (numeric matrix) with dimensions \code{d} by \code{r}.
#' @param dep The true dependence parameter (numeric scalar).
#'
#' @return A numeric scalar representing the standardized mean squared error combining the matrix and dependence parameter estimation errors.
#'
#' @details
#' The function first ensures each estimated matrix has the same number of columns as the true matrix by padding with zeros if necessary.
#' It then calculates the average normalized squared differences between the estimated and true matrices, scaling by the maximum estimated values across seeds to avoid division by zero.
#' The dependence parameter error is incorporated as a standardized squared difference.
#'
#' @references
#' Mourahib, A., Kiriliouk, A., & Segers, J. (2025). A penalized least squares estimator for extreme-value mixture models. \emph{arXiv preprint} arXiv:2506.15272.
#'
#' @export

SMSE_log <- function (list, A , dep){
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



#' Extract Upper Triangular Values from a Matrix as a Vector
#'
#' Extracts the elements above the main diagonal (upper triangular part) of a square matrix as a vector.
#' Optionally, specified indices can be excluded from the extracted values.
#'
#' @param mat A numeric square matrix.
#' @param ind Optional integer vector of indices to exclude from the extracted upper triangular values.
#'
#' @return A numeric vector containing the upper triangular elements of \code{mat}, with specified indices removed if provided.
#'
#' @export


matrix_vector <- function(mat, ind = NULL){
  upper_values <- mat[upper.tri(mat, diag = FALSE)]
  if(!is.null(ind)){
    upper_values <- upper_values[-ind]
  }
  return(upper_values)
}

#' Compute standardized Mean Squared Error in Equation (4.4) from Mourahib, Kiriliouk, and Segers (2025) for the mixture Hüsler–Reiss model
#'
#' Calculates a standardized mean squared error (SMSE) between a list of estimated matrices and a true coefficient matrix,
#' incorporating both matrix and dependence parameter errors. Scaling adjusts for variability across estimates.
#'
#' @param list_matrices A list where each element contains an \code{Estimation} list with \code{pls_matrix} and \code{pls_dep} components.
#' @param A The true coefficient matrix.
#' @param ind Optional integer vector of indices to exclude when comparing dependence parameters.
#'
#' @return A numeric value representing the aggregated standardized mean squared error over all estimates.
#'
#' @references
#' Mourahib, A., Kiriliouk, A., & Segers, J. (2025). A penalized least squares estimator for extreme-value mixture models. \emph{arXiv preprint} arXiv:2506.15272.
#'
#' @export

SMSE_HR <- function (list_matrices , A , ind = NULL  ){
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

