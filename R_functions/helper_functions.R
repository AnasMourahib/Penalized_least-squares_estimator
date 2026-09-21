extract_signatures <- function(A) {
  if (!is.matrix(A)) stop("A must be a matrix.")

  signatures <- lapply(seq_len(ncol(A)), function(k) which(A[, k] > 0))
  unique(signatures)
}


EDS <- function(L1, L2) {

  s1 <- sapply(L1, function(x) paste(sort(x), collapse = ","))
  s2 <- sapply(L2, function(x) paste(sort(x), collapse = ","))

  inters <- intersect(s1, s2)
  uni <- union(s1, s2)

  1 - length(inters) / length(uni)
}



# Puts the columns of start in the "right" order (i.e. the order corresponding to A)
# works also if ncol(start) > ncol(A)
simulate_sparse_A <- function(
    d,
    r,
    sparse_proportion,
    min_nonzero_per_row = 2L,
    require_nonzero_columns = TRUE,
    max_attempts = 10000L
) {

  # ------------------------------------------------------------
  # 1. Validate the inputs
  # ------------------------------------------------------------

  if (d < 1L || r < 1L) {
    stop("d and r must be positive integers.")
  }

  if (
    sparse_proportion < 0 ||
    sparse_proportion >= 1
  ) {
    stop("sparse_proportion must be in [0, 1).")
  }

  if (
    min_nonzero_per_row < 1L ||
    min_nonzero_per_row > r
  ) {
    stop(
      "min_nonzero_per_row must be between 1 and r."
    )
  }

  # If every column must be nonzero, at most 2^d - 1
  # distinct nonempty support patterns are available.
  if (
    require_nonzero_columns &&
    r > 2^d - 1L
  ) {
    stop(
      "Cannot generate ",
      r,
      " distinct nonempty column signatures with d = ",
      d,
      ". The maximum is 2^d - 1 = ",
      2^d - 1L,
      "."
    )
  }

  # ------------------------------------------------------------
  # 2. Determine the requested number of nonzero entries
  # ------------------------------------------------------------

  total_entries <- d * r

  requested_zeros <- round(
    sparse_proportion * total_entries
  )

  target_nonzero <- total_entries - requested_zeros

  # Row constraint requires at least this many nonzero entries.
  minimum_from_rows <- d * min_nonzero_per_row

  # Nonzero-column constraint requires at least r entries.
  minimum_from_columns <- if (require_nonzero_columns) {
    r
  } else {
    0L
  }

  minimum_nonzero <- max(
    minimum_from_rows,
    minimum_from_columns
  )

  if (target_nonzero < minimum_nonzero) {
    stop(
      "The requested sparsity is incompatible with the constraints. ",
      "The requested number of nonzero entries is ",
      target_nonzero,
      ", but at least ",
      minimum_nonzero,
      " entries are required."
    )
  }

  # ------------------------------------------------------------
  # 3. Generate a binary support matrix
  # ------------------------------------------------------------

  valid_support <- FALSE
  support_matrix <- NULL

  for (attempt in seq_len(max_attempts)) {

    support_candidate <- matrix(
      0L,
      nrow = d,
      ncol = r
    )

    selected_positions <- sample.int(
      total_entries,
      size = target_nonzero,
      replace = FALSE
    )

    support_candidate[selected_positions] <- 1L

    # Every row must contain at least the requested number
    # of nonzero entries.
    valid_rows <- all(
      rowSums(support_candidate) >= min_nonzero_per_row
    )

    # Optionally ensure every direction is nonempty.
    valid_columns <- (
      !require_nonzero_columns ||
        all(colSums(support_candidate) > 0L)
    )

    # Convert each column support into a signature.
    column_signatures <- apply(
      support_candidate,
      2,
      paste0,
      collapse = ""
    )

    unique_columns <- (
      length(unique(column_signatures)) == r
    )

    if (
      valid_rows &&
      valid_columns &&
      unique_columns
    ) {
      support_matrix <- support_candidate
      valid_support <- TRUE
      break
    }
  }

  if (!valid_support) {
    stop(
      "A valid support matrix was not found after ",
      max_attempts,
      " attempts. The combination of d, r, sparsity, ",
      "and row constraints may be difficult or infeasible."
    )
  }

  # ------------------------------------------------------------
  # 4. Generate positive values on the selected support
  # ------------------------------------------------------------

  A <- matrix(
    0,
    nrow = d,
    ncol = r
  )

  number_nonzero <- sum(support_matrix)

  A[support_matrix == 1L] <- runif(
    number_nonzero,
    min = 0.1,
    max = 0.9
  )

  # ------------------------------------------------------------
  # 5. Normalize every row to sum to one
  # ------------------------------------------------------------

  A <- A / rowSums(A)

  # ------------------------------------------------------------
  # 6. Defensive checks
  # ------------------------------------------------------------

  estimated_signatures <- apply(
    A > 0,
    2,
    paste0,
    collapse = ""
  )

  if (length(unique(estimated_signatures)) != r) {
    stop("Internal error: duplicated column signatures.")
  }

  if (
    any(
      rowSums(A > 0) < min_nonzero_per_row
    )
  ) {
    stop("Internal error: a row has too few nonzero entries.")
  }

  if (!isTRUE(all.equal(
    rowSums(A),
    rep(1, d),
    tolerance = 1e-12
  ))) {
    stop("Internal error: rows do not sum to one.")
  }

  return(A)
}



empirical_cdf <- function(x) {
  # rank each x[i] by how many values are ≤ it,
  # then divide by n to get the empirical CDF at x[i]
  rank(x, ties.method = "max") / length(x)
}



fun_estimate_empirical_corr <- function(X , q){
  N <- nrow(X)
  matrix_chi_emp <- matrix( 0 , nrow = d , ncol = d)
  for( s in 1:d){
    vec_s <- X[,s]
    for(t in 1:d){
      vec_t <- X[,t]
      cdf_s <- empirical_cdf(vec_s)
      cdf_t <- empirical_cdf(vec_t)
      chi_st <- (length(which( (cdf_s >q)  &    (cdf_t >q) )    ) / N )   / (1 - q)
      if(chi_st >1) {chi_st <- 1}
      matrix_chi_emp[s,t] <- chi_st
    }
  }
  chi_emp <- matrix_chi_emp[upper.tri(matrix_chi_emp)]
  return(chi_emp)
}


Noise_simulator <- function( N , d , shape) {
  bool = FALSE
  while(bool == FALSE) {
    df = (d^2 - d)/2
    Theta = runif(df , 0 , 0.1)
    R = p2P(Theta)
    if(is.positive.definite(R)) {
      bool = TRUE
    }
  }
  Theta = P2p(R)
  #print(Theta)
  myCop <- normalCopula(param = Theta, dim = d, dispstr = "un")
  myMvd <- mvdc(
    copula = myCop,
    margins = rep("unif", d),
    paramMargins = replicate(d, list(min = 0, max = 1), simplify = FALSE)
  )
  U <- rMvdc(N, myMvd)
  frechet_quantile <- function(u) {
    (-log(u))^(-1 / shape)
  }

  X <- apply(U , c(1 , 2) , frechet_quantile)
  return(X)
}

Noise_simulator_ind <- function(N, d, shape) {

  # Independent uniforms
  U <- matrix(runif(N * d), nrow = N, ncol = d)

  # Frechet quantile transformation
  X <- (-log(U))^(-1 / shape)

  return(X)
}


shuffleCols <- function(start, A){
  r <- ncol(A)
  k <- ncol(start)
  perms <- permutations(n = k, r = r)
  temp <- apply(perms, 1, function(j) sum(abs(start[,j] - A)))
  indx <- which(temp == min(temp))
  firstcols <- perms[indx,]
  lastcols <- setdiff(c(1:k),firstcols)
  return(start[,c(firstcols,lastcols)])
}


construct_symmetric_matrix_2 <- function(vec) {
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





cross_validation<-function(d , r, A , grid, lambda, num_col = NULL, start , type = c("SSR_row_HR", "SSR_row_log"), p , w , num_class=10){
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
                   as.integer(num_col), as.integer(q), rep(theta_alpha, num_col), as.double(w),
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
  start_dep <- if (v == 1) runif(1, 0.1, 0.9) else rep(1,v)
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

  # Partial Estimation of ?_Z
  interm_3 <- function(theta) {
    if (any(theta < 0) || (type == "SSR_row_log" && any(theta > 1))) return(1e16)
    result <- .C(type, as.double(p), as.double(lambda), A_vec_fixed, as.integer(d),
                 as.integer(num_col), as.integer(q), as.double(theta), as.double(w),
                 grid_flat, R = double(1))$R
    if (!is.finite(result)) return(1e16)
    result
  }

  temp <- optim(theta_Z_pr, interm_3, method = "L-BFGS-B", lower = rep(0, v), control = list(maxit = 1000))
  theta_Z_final <- if (type == "SSR_row_HR") construct_symmetric_matrix_2(temp$par) else temp$par

  return(list(pls_matrix = matrix_A_PLS, pls_dep = theta_Z_final))
}





param_estim_2 <- function(d , r  , grid , lambda , num_col = NULL , start , type = c("SSR_row_HR", "SSR_row_log"), p ,  w , task = NULL , seed ){
  set.seed(seed)
  q <- nrow(grid)
  if(is.null(num_col)){num_col <- r} #if not specified, use the correct number of columns
  l <- d*num_col
  v <- 1
  interm <- function(theta){
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
                    as.integer(num_col), as.integer(q),  as.double(rep(theta_Gamma , (d^(2) -d )/2  )), as.double(w), as.double(c(t(grid))), R = double(1))$R
    }
    if (!is.finite(result)) {
      return(10^16)  # Return a large penalty value to handle non-finite values
    }
    return(result)
  }
  if(v == 1){
    start_dep <- runif(1 , min=0.3 , max = 0.6)
  }
  else{
    start_dep <- runif( v  , min = 0.2 , max = 2)
  }
  start <- c(start , start_dep)
  temp <- optim( start , interm , method="L-BFGS-B" , lower = c( rep(0, l  ) , rep(0.3  , v) ) , upper =  c( rep(Inf, l  ) , rep(0.6  , v) )  , control = list('maxit' = 1000))
  estim_pr <- temp$par
  theta_A_pr <- estim_pr[1:l]  #This is the estimation asscociated to the matrix
  A_pr <- matrix(   theta_A_pr , nrow = d , byrow = TRUE  )
  theta_Z_pr <- estim_pr[(l+1) : (l+v)]
  if(type == "SSR_row_HR"){
    theta_Z_pr <- construct_symmetric_matrix_2(rep(temp$par , (d^2 - d)/2  ))
  }

  #####If we just need to estimate the extreme directions then we can stop here as we do not want to compute SMSE
  if(task=="ED_identification"){
    return(list("pls_matrix" =  A_pr ,  "pls_dep" = theta_Z_pr ))
  }

  else{
    ################ Fix $\Theta_Z$ and estimate A partially

    interm_2 <- function(theta ){
      if( any(  theta >1 | theta <0   ) ){
        return(10^16)
      }
      if (type == "SSR_row_log"){
        result <- .C(type , as.double(p),as.double(lambda), as.double(theta), as.integer(d),
                     as.integer(num_col), as.integer(q), as.double(rep(theta_Z_pr, num_col)), as.double(w), as.double(c(t(grid))), R = double(1))$R
      }
      if(type == "SSR_row_HR"){
        result <- .C(type , as.double(p),as.double(lambda), as.double(theta), as.integer(d),
                     as.integer(num_col), as.integer(q), as.double(rep(theta_Z_pr , (d^(2) -d )/2  )), as.double(w), as.double(c(t(grid))), R = double(1))$R
      }
      if (!is.finite(result)) {
        return(10^16)  # Return a large penalty value to handle non-finite values
      }
      return(result)
    }

    start <- theta_A_pr
    temp <- optim( start , interm_2  , method="L-BFGS-B" , lower=rep(0 , l) , control = list('maxit' = 1000) )
    PLS_A <- temp$par
    matrix_A_PLS <- matrix( normalize_group( PLS_A , num_col )  , ncol = num_col , byrow = T)

    ############## Fix A and estimate $\Theta_Z$ partially
    interm_3 <- function(theta ){
      if( any( theta <0   ) ){
        return(10^16)
      }
      if (type == "SSR_row_log"){
        if( any( theta > 1   ) ){
          return(10^16)
        }
        ######At this stage, since $A$ is fixed, you can choose lambda = 0
        result <- .C(type , as.double(p),as.double(lambda), as.double(as.vector(t(matrix_A_PLS))), as.integer(d),
                     as.integer(num_col), as.integer(q), as.double(rep(theta, num_col)), as.double(w), as.double(c(t(grid))), R = double(1))$R
      }
      if(type == "SSR_row_HR"){
        result <- .C(type , as.double(p),as.double(lambda), as.double(as.vector(t(matrix_A_PLS))), as.integer(d),
                     as.integer(num_col), as.integer(q), as.double(rep(theta , (d^(2) -d )/2  )), as.double(w), as.double(c(t(grid))), R = double(1))$R
      }
      if (!is.finite(result)) {
        return(10^16)  # Return a large penalty value to handle non-finite values
      }
      return(result)
    }
    start <- theta_Z_pr
    temp <- optim( start , interm_3 , method="L-BFGS-B" , lower=rep(0 , v) , control = list('maxit' = 1000) )
    if(type == "SSR_row_HR"){
      theta_Z_pr <- construct_symmetric_matrix_2(rep(temp$par , (d^2 - d)/2  ))
    }
    return(list("pls_matrix" =  matrix_A_PLS ,  "pls_dep" = theta_Z_pr ))
  }
}



#####This param_estim and cross validation to parallelize the folds so that we use a warm start instead of paralellizing over the lambda's


param_estim_path_fold <- function(
    d,
    r,
    grid,
    lambda,
    num_col = NULL,
    start,
    start_dep = NULL,
    type = c("SSR_row_HR", "SSR_row_log"),
    p,
    w,
    task = NULL,
    seed = NULL,
    maxit = 1000
) {

  type <- match.arg(type)

  if (is.null(num_col)) {
    num_col <- r
  }

  q <- nrow(grid)
  l <- d * num_col
  v <- 1L

  # Deterministic initial value for the dependence parameter
  if (is.null(start_dep)) {
    start_dep <- runif(1 , 0.15 , 0.25)
  }

  # Make sure the starting values satisfy the bounds
  start <- pmin(pmax(start, 0), 1)
  #start_dep <- pmin(pmax(start_dep, 0.3), 0.6)

  initial_par <- c(start, start_dep)

  # Precompute constant C inputs
  p_C <- as.double(p)
  lambda_C <- as.double(lambda)
  d_C <- as.integer(d)
  num_col_C <- as.integer(num_col)
  q_C <- as.integer(q)
  w_C <- as.double(w)
  grid_C <- as.double(t(grid))

  interm <- function(theta) {

    theta_A <- theta[seq_len(l)]
    theta_dep <- theta[l + 1L]

    if (type == "SSR_row_log") {

      dependence_C <- as.double(
        rep(theta_dep, num_col)
      )

    } else {

      dependence_C <- as.double(
        rep(theta_dep, (d^2 - d) / 2)
      )
    }

    result <- .C(
      type,
      p_C,
      lambda_C,
      as.double(theta_A),
      d_C,
      num_col_C,
      q_C,
      dependence_C,
      w_C,
      grid_C,
      R = double(1)
    )$R

    if (!is.finite(result)) {
      return(1e16)
    }

    result
  }

  temp <- optim(
    par = initial_par,
    fn = interm,
    method = "L-BFGS-B",
    lower = c(
      rep(0, l),
      rep(0.1, v)
    ),
    upper = c(
      rep(1, l),
      rep(0.3, v)
    ),
    control = list(
      maxit = maxit
    )
  )

  estimated_parameters <- temp$par

  theta_A_pr <- estimated_parameters[seq_len(l)]

  A_pr <- matrix(
    theta_A_pr,
    nrow = d,
    ncol = num_col,
    byrow = TRUE
  )

  theta_dep_pr <- estimated_parameters[l + 1L]

  if (type == "SSR_row_HR") {

    theta_Z_pr <- construct_symmetric_matrix_2(
      rep(theta_dep_pr, (d^2 - d) / 2)
    )

  } else {

    theta_Z_pr <- theta_dep_pr
  }

  # Current task stops after the first optimization
  if (identical(task, "ED_identification")) {

    return(list(
      pls_matrix = A_pr,
      pls_dep = theta_Z_pr,

      # Additional outputs required for warm starts
      par = temp$par,
      convergence = temp$convergence,
      objective = temp$value,
      counts = temp$counts,
      message = temp$message
    ))
  }

  # Keep the remaining part of your original function here
  # if you use tasks other than ED_identification.

  stop(
    "The optimized version currently covers task = 'ED_identification'. ",
    "Add your original second and third optimization stages here ",
    "for other tasks."
  )
}



cross_validation_path_fold <- function(
    class_k,
    lambda_grid,
    d,
    r,
    grid,
    num_col = NULL,
    start,
    type = c("SSR_row_HR", "SSR_row_log"),
    p,
    w_train,
    w_test,
    task,
    seed,
    maxit_cv = 250
) {

  type <- match.arg(type)

  if (is.null(num_col)) {
    num_col <- r
  }

  q <- nrow(grid)
  l <- d * num_col

  # Start with the largest penalty
  lambda_order <- order(
    lambda_grid,
    decreasing = TRUE
  )

  lambda_sorted <- lambda_grid[lambda_order]

  scores_sorted <- rep(
    Inf,
    length(lambda_sorted)
  )

  # Initial values for the first lambda
  current_start_A <- start
  current_start_dep <- 0.45

  # Precompute constant C arguments
  p_C <- as.double(p)
  zero_C <- as.double(0)
  d_C <- as.integer(d)
  num_col_C <- as.integer(num_col)
  q_C <- as.integer(q)
  grid_C <- as.double(t(grid))
  test_w_C <- as.double(w_test[[class_k]])

  for (j in seq_along(lambda_sorted)) {

    lambda_j <- lambda_sorted[j]

    fit <- tryCatch(
      param_estim_path_fold(
        d = d,
        r = r,
        grid = grid,
        lambda = lambda_j,
        num_col = num_col,
        start = current_start_A,
        start_dep = current_start_dep,
        type = type,
        p = p,
        w = w_train[[class_k]],
        task = task,
        seed = seed,
        maxit = maxit_cv
      ),
      error = function(e) {
        NULL
      }
    )

    if (is.null(fit)) {
      next
    }

    v_A <- as.double(
      as.vector(t(fit$pls_matrix))
    )

    if (type == "SSR_row_log") {

      dependence_C <- as.double(
        rep(fit$pls_dep, num_col)
      )

    } else {

      # fit$pls_dep is already the dependence structure
      dependence_C <- as.double(fit$pls_dep)
    }

    test_score <- .C(
      type,
      p_C,
      zero_C,
      v_A,
      d_C,
      num_col_C,
      q_C,
      dependence_C,
      test_w_C,
      grid_C,
      R = double(1)
    )$R

    if (is.finite(test_score)) {
      scores_sorted[j] <- test_score
    }

    # Warm start only if optimization produced valid parameters
    if (
      !is.null(fit$par) &&
      length(fit$par) == l + 1L &&
      all(is.finite(fit$par))
    ) {

      current_start_A <- fit$par[seq_len(l)]
      current_start_dep <- fit$par[l + 1L]
    }
  }

  # Restore the original lambda ordering
  scores <- rep(
    Inf,
    length(lambda_grid)
  )

  scores[lambda_order] <- scores_sorted

  scores
}


########



cross_validation_2<-function(d , r , grid, lambda, num_col = NULL, start , type = c("SSR_row_HR", "SSR_row_log"), p , w , num_class=5, task , seed){
  #start <- c(start, 0.5)
  w_train <- w$train
  w_test <- w$test
  q <- nrow(grid)
  if(is.null(num_col)){num_col <- r} #if not specified, use the correct number of columns
  l <- d * num_col
  CV <- vector(length = num_class)
  for (class_k in 1:num_class){
    optimizer_minus_class_k <- param_estim_2(d = d , r = r  , grid = grid , lambda = lambda , num_col = num_col , start = start  , type = type , p = p  ,  w = w_train[[class_k]] , task , seed = seed)

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









N_generate_Mix_log <- function(N, A, alpha) {
  r <- ncol(A)
  if(length(alpha) == 1){alpha <- rep(alpha,r)}
  Z <- vector('list', length = r)
  for(k  in 1:r){
    Z[[k]] <- matrix(0,nrow(A),N)
    sig <- which(A[,k]>0)
    Z[[k]][sig,] <- t(rmev(N, length(sig), param = 1/alpha[k], model = "log"))*A[sig,k]
  }
  M <- apply(simplify2array(Z),c(1,2),max)
  return(t(M))
}



# ? = sqrt(G)/2.
# The relationship between ? and r�the dependence parameter of the H�sler�Reiss (HR) model�
# in the R package `evd` is discussed in the documentation of the `rmev` function.
# Express r in terms of G using the stdf expression
# provided in Section 4.2.2 of https://link.springer.com/article/10.1007/s10687-024-00501-4
# and the details given on page 14 of http://cran.fhcrc.org/web/packages/evd/evd.pdf.
N_generate_Mix_hr <- function(N, A, sigma) {
  r <- ncol(A)
  Z <- vector('list', length = r)
  for(k  in 1:r){
    Z[[k]] <- matrix(0,nrow(A),N)
    sig <- which(A[,k]>0)
    lsig <- length(sig)
    if(lsig == 1){ Z[[k]][sig , ] <- -1/ log( runif(N , 0 , 1) ) * A[sig,k]  }
    else{
      sub_sigma <- sigma[sig , sig]
      Z[[k]][sig,] <- t(rmev(N, lsig, sigma = sub_sigma, model = "hr"))*A[sig,k]
    }
  }
  M <- apply(simplify2array(Z),c(1,2),max)
  return(t(M))
}







starting_point <-function(data, nrcol, quant = 0.9){
  N <- nrow(data)

  dataP <- apply(data, 2, function(i) N/(N + 0.5 - rank(i)))
  U <- quantile(rowSums(dataP),quant)
  dataU <-dataP[rowSums(dataP)>U,]
  ndata <- t(apply(dataU,1, function(i) i/sum(i)))

  kmean <- kmeans(ndata,centers=nrcol,nstart=5)
  startk <- sapply(c(1:nrcol), function(j) kmean$centers[j,]*kmean$size[j])
  resk <- t(apply(startk, 1, function(x) x/sum(x)))
  return(resk)
}
