param_estim <- function(d, r, A, grid, lambda, num_col = NULL, start, type = c("SSR_row_HR", "SSR_row_log"), p, w) {
  q <- nrow(grid)
  if (is.null(num_col)) num_col <- r
  l <- d * num_col
  v <- ifelse(type == "SSR_row_log", 1, (d * (d - 1)) / 2)
  
  grid_flat <- as.double(c(t(grid)))  # precompute
  A_length <- d * num_col
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
  A_vec_fixed <- as.vector(t(matrix_A_PLS))
  
  # Partial Estimation of θ_Z
  interm_3 <- function(theta) {
    if (any(theta < 0) || (type == "SSR_row_log" && any(theta > 1))) return(1e16)
    result <- .C(type, as.double(p), as.double(lambda), A_vec_fixed, as.integer(d),
                 as.integer(num_col), as.integer(q), as.double(theta), as.double(w),
                 grid_flat, R = double(1))$R
    if (!is.finite(result)) return(1e16)
    result
  }
  
  temp <- optim(theta_Z_pr, interm_3, method = "L-BFGS-B", lower = rep(0, v), control = list(maxit = 1000))
  theta_Z_final <- if (type == "SSR_row_HR") construct_symmetric_matrix(temp$par) else temp$par
  
  return(list(pls_matrix = matrix_A_PLS, pls_dep = theta_Z_final))
}
