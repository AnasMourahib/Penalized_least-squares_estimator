# ============================================================
# Helper functions required by estimation_main.R
# stdfEmp() must be supplied by its source file or package.
# ============================================================

construct_symmetric_matrix_2 <- function(vec) {
  # Determine d from the length of vec
  d <- (1 + sqrt(1 + 8 * length(vec))) / 2
  if (d != floor(d)) stop("Vector length is incorrect for a symmetric matrix")
  d <- as.integer(d)

  # Initialize d x d matrix with 1s on the diagonal
  mat <- diag(  0, d, d)

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
  if (length(v) %% group_size != 0L) stop("length(v) must be divisible by group_size.")
  out <- v
  for (i in seq(1L, length(v), by = group_size)) {
    index <- i:(i + group_size - 1L); total <- sum(v[index])
    if (total > 0) out[index] <- v[index] / total
  }
  out
}

W_calculus <- function(k, num_class, X, grid, q = nrow(grid)) {
  if (q != nrow(grid)) stop("q must equal nrow(grid).")
  N <- nrow(X); folds <- split(seq_len(N), cut(seq_len(N), num_class, labels = FALSE))
  W_train <- W_test <- vector("list", num_class)
  for (class_k in seq_len(num_class)) {
    test_index <- folds[[class_k]]; train_index <- setdiff(seq_len(N), test_index)
    R_train <- apply(X[train_index, , drop = FALSE], 2, rank)
    R_test <- apply(X[test_index, , drop = FALSE], 2, rank)
    k_train <- max(1, round(length(train_index) / N * k)); k_test <- max(1, round(length(test_index) / N * k))
    W_train[[class_k]] <- vapply(seq_len(q), function(m) stdfEmp(R_train, k_train, grid[m, ]), numeric(1))
    W_test[[class_k]] <- vapply(seq_len(q), function(m) stdfEmp(R_test, k_test, grid[m, ]), numeric(1))
  }
  list(train = W_train, test = W_test)
}

starting_point <- function(data, nrcol, quant = 0.9) {
  N <- nrow(data); dataP <- apply(data, 2, function(x) N / (N + 0.5 - rank(x)))
  threshold <- stats::quantile(rowSums(dataP), quant); dataU <- dataP[rowSums(dataP) > threshold, , drop = FALSE]
  angular <- t(apply(dataU, 1, function(x) x / sum(x)))
  if (nrow(angular) < nrcol) stop("Not enough threshold exceedances for the requested number of columns.")
  fit <- stats::kmeans(angular, centers = nrcol, nstart = 5)
  centers <- sapply(seq_len(nrcol), function(j) fit$centers[j, ] * fit$size[j])
  t(apply(centers, 1, function(x) x / sum(x)))
}

extract_signatures <- function(A) {
  if (!is.matrix(A)) stop("A must be a matrix.")
  unique(lapply(seq_len(ncol(A)), function(k) which(A[, k] > 0)))
}

EDS <- function(L1, L2) {
  s1 <- vapply(L1, function(x) paste(sort(x), collapse = ","), character(1))
  s2 <- vapply(L2, function(x) paste(sort(x), collapse = ","), character(1))
  union_set <- union(s1, s2); if (!length(union_set)) return(0)
  1 - length(intersect(s1, s2)) / length(union_set)
}

simulate_sparse_A <- function(d, r, sparse_proportion, min_nonzero_per_row = 2L,
                              require_nonzero_columns = TRUE, max_attempts = 10000L) {
  if (d < 1L || r < 1L || sparse_proportion < 0 || sparse_proportion >= 1) stop("Invalid dimensions or sparsity.")
  target <- d * r - round(sparse_proportion * d * r)
  minimum <- max(d * min_nonzero_per_row, if (require_nonzero_columns) r else 0L)
  if (target < minimum || (require_nonzero_columns && r > 2^d - 1L)) stop("Requested support is infeasible.")
  for (attempt in seq_len(max_attempts)) {
    S <- matrix(0L, d, r); S[sample.int(d * r, target)] <- 1L
    signatures <- apply(S, 2, paste0, collapse = "")
    valid <- all(rowSums(S) >= min_nonzero_per_row) && (!require_nonzero_columns || all(colSums(S) > 0L)) && length(unique(signatures)) == r
    if (valid) {
      A <- matrix(0, d, r); A[S == 1L] <- runif(sum(S), 0.1, 0.9); return(A / rowSums(A))
    }
  }
  stop("No valid support found.")
}

empirical_cdf <- function(x) rank(x, ties.method = "max") / length(x)

fun_estimate_empirical_corr <- function(X, q) {
  d <- ncol(X); N <- nrow(X); U <- apply(X, 2, empirical_cdf); chi <- matrix(0, d, d)
  for (s in seq_len(d)) for (t in seq_len(d)) chi[s, t] <- min(1, mean(U[, s] > q & U[, t] > q) / (1 - q))
  chi[upper.tri(chi)]
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
N_generate_Mix_hr <- function(N, A, Gamma) {
  sigma <- Gamma /4
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

