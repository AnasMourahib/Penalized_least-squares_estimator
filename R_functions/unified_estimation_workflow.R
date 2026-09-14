# ============================================================
# 1. Unified parameter estimation
# ============================================================
param_estim <- function(d, r, grid, lambda, num_col = NULL, start, start_dep = NULL,
                        type = c("SSR_row_HR", "SSR_row_log"), p, w, task = NULL,
                        seed = NULL, maxit = 1000) {
  type <- match.arg(type)
  if (is.null(num_col)) num_col <- r
  if (!is.null(seed)) set.seed(seed)
  q <- nrow(grid); l <- d * num_col; v <- 1L
  if (length(start) != l) stop("start must have length d * num_col.")
  if (is.null(start_dep)) start_dep <- runif(1, 0.15, 0.25)
  start <- pmin(pmax(start, 0), 1); start_dep <- pmin(pmax(start_dep, 0.1), 0.3)

  p_C <- as.double(p); lambda_C <- as.double(lambda); d_C <- as.integer(d)
  num_col_C <- as.integer(num_col); q_C <- as.integer(q)
  w_C <- as.double(w); grid_C <- as.double(t(grid))

  objective_joint <- function(theta) {
    theta_A <- theta[seq_len(l)]; theta_dep <- theta[l + 1L]
    dep_C <- if (type == "SSR_row_log") as.double(rep(theta_dep, num_col)) else as.double(rep(theta_dep, (d^2 - d) / 2))
    value <- .C(type, p_C, lambda_C, as.double(theta_A), d_C, num_col_C, q_C, dep_C, w_C, grid_C, R = double(1))$R
    if (!is.finite(value)) 1e16 else value
  }

  fit_joint <- optim(c(start, start_dep), objective_joint, method = "L-BFGS-B",
                     lower = c(rep(0, l), rep(0.1, v)), upper = c(rep(1, l), rep(0.3, v)),
                     control = list(maxit = maxit))
  theta_A <- fit_joint$par[seq_len(l)]; theta_dep <- fit_joint$par[l + 1L]
  A_joint <- matrix(theta_A, nrow = d, ncol = num_col, byrow = TRUE)
  dep_joint <- if (type == "SSR_row_HR") construct_symmetric_matrix_2(rep(theta_dep, (d^2 - d) / 2)) else theta_dep

  if (identical(task, "ED_identification")) {
    return(list(pls_matrix = A_joint, pls_dep = dep_joint, par = fit_joint$par,
                convergence = fit_joint$convergence, objective = fit_joint$value,
                counts = fit_joint$counts, message = fit_joint$message))
  }

  objective_A <- function(theta) {
    if (any(theta > 1 | theta < 0)) return(1e16)
    dep_C <- if (type == "SSR_row_log") as.double(rep(theta_dep, num_col)) else as.double(rep(theta_dep, (d^2 - d) / 2))
    value <- .C(type, p_C, lambda_C, as.double(theta), d_C, num_col_C, q_C, dep_C, w_C, grid_C, R = double(1))$R
    if (!is.finite(value)) 1e16 else value
  }

  fit_A <- optim(theta_A, objective_A, method = "L-BFGS-B", lower = rep(0, l), control = list(maxit = maxit))
  A_final <- matrix(normalize_group(fit_A$par, num_col), ncol = num_col, byrow = TRUE)

  objective_dep <- function(theta) {
    if (any(theta < 0) || (type == "SSR_row_log" && any(theta > 1))) return(1e16)
    dep_C <- if (type == "SSR_row_log") as.double(rep(theta, num_col)) else as.double(rep(theta, (d^2 - d) / 2))
    value <- .C(type, p_C, lambda_C, as.double(as.vector(t(A_final))), d_C, num_col_C, q_C, dep_C, w_C, grid_C, R = double(1))$R
    if (!is.finite(value)) 1e16 else value
  }

  fit_dep <- optim(theta_dep, objective_dep, method = "L-BFGS-B", lower = rep(0, v), control = list(maxit = maxit))
  dep_final <- if (type == "SSR_row_HR") construct_symmetric_matrix_2(rep(fit_dep$par, (d^2 - d) / 2)) else fit_dep$par
  list(pls_matrix = A_final, pls_dep = dep_final, par = c(as.vector(t(A_final)), fit_dep$par),
       convergence = fit_dep$convergence, objective = fit_dep$value, counts = fit_dep$counts,
       message = fit_dep$message, joint_fit = fit_joint, A_fit = fit_A, dependence_fit = fit_dep)
}

# ============================================================
# 2A. ED identification: one fold, complete lambda path, warm starts
# ============================================================
cross_validation_path_fold <- function(class_k, lambda_grid, d, r, grid, num_col = NULL,
                                       start, type = c("SSR_row_HR", "SSR_row_log"), p,
                                       w_train, w_test, task, seed, maxit_cv = 250) {
  type <- match.arg(type)
  if (is.null(num_col)) num_col <- r
  q <- nrow(grid); l <- d * num_col
  lambda_order <- order(lambda_grid, decreasing = TRUE); lambda_sorted <- lambda_grid[lambda_order]
  scores_sorted <- rep(Inf, length(lambda_sorted)); current_start_A <- start; current_start_dep <- 0.2
  p_C <- as.double(p); zero_C <- as.double(0); d_C <- as.integer(d)
  num_col_C <- as.integer(num_col); q_C <- as.integer(q); grid_C <- as.double(t(grid))
  test_w_C <- as.double(w_test[[class_k]])

  for (j in seq_along(lambda_sorted)) {
    fit <- tryCatch(param_estim(d, r, grid, lambda_sorted[j], num_col, current_start_A,
                              current_start_dep, type, p, w_train[[class_k]], task, seed, maxit_cv),
                    error = function(e) NULL)
    if (is.null(fit)) next
    v_A <- as.double(as.vector(t(fit$pls_matrix)))
    dep_C <- if (type == "SSR_row_log") as.double(rep(fit$pls_dep, num_col)) else as.double(fit$pls_dep)
    score <- .C(type, p_C, zero_C, v_A, d_C, num_col_C, q_C, dep_C, test_w_C, grid_C, R = double(1))$R
    if (is.finite(score)) scores_sorted[j] <- score
    if (!is.null(fit$par) && length(fit$par) == l + 1L && all(is.finite(fit$par))) {
      current_start_A <- fit$par[seq_len(l)]; current_start_dep <- fit$par[l + 1L]
    }
  }
  scores <- rep(Inf, length(lambda_grid)); scores[lambda_order] <- scores_sorted; scores
}

# ============================================================
# 2B. Other tasks: one lambda, all folds, no warm starts
# ============================================================
cross_validation_standard <- function(lambda, d, r, grid, num_col = NULL, start,
                                      type = c("SSR_row_HR", "SSR_row_log"), p, w,
                                      num_class = 5, task, seed, maxit_cv = 250) {
  type <- match.arg(type)
  if (is.null(num_col)) num_col <- r
  q <- nrow(grid); scores <- numeric(num_class)
  for (class_k in seq_len(num_class)) {
    fit <- param_estim(d, r, grid, lambda, num_col, start, NULL, type, p,
                       w$train[[class_k]], task, seed, maxit_cv)
    v_A <- as.double(as.vector(t(fit$pls_matrix)))
    dep_C <- if (type == "SSR_row_log") as.double(rep(fit$pls_dep, num_col)) else as.double(fit$pls_dep)
    scores[class_k] <- .C(type, as.double(p), as.double(0), v_A, as.integer(d),
                          as.integer(num_col), as.integer(q), dep_C,
                          as.double(w$test[[class_k]]), as.double(t(grid)), R = double(1))$R
  }
  mean(scores)
}

# ============================================================
# 3. Main dispatcher
# ED: parallel folds + warm starts. Other tasks: parallel lambdas.
# refined_grid = 0 means broad grid only; otherwise use refinement.
# ============================================================
main_fit <- function(X, w, w_total, lambda_grid, grid, num_col = NULL, start = NULL,
                     type = c("SSR_row_HR", "SSR_row_log"), k, p, num_class = 5, cl,
                     d, r, task, seed, maxit_cv = 250, maxit_final = 1000,
                     cv_tolerance = 0.001, refined_grid_length = 15, refined_grid = 0) {
  type <- match.arg(type)
  if (is.null(num_col)) num_col <- r
  if (is.null(start)) start <- as.vector(t(starting_point(X, num_col)))
  if (length(start) != d * num_col) stop("start must have length d * num_col.")

  evaluate_grid <- function(current_grid) {
    if (identical(task, "ED_identification")) {
      fold_scores <- parallel::parLapply(cl, seq_len(num_class), cross_validation_path_fold,
        lambda_grid = current_grid, d = d, r = r, grid = grid, num_col = num_col,
        start = start, type = type, p = p, w_train = w$train, w_test = w$test,
        task = task, seed = seed, maxit_cv = maxit_cv)
      score_matrix <- do.call(rbind, fold_scores)
      scores <- colMeans(score_matrix)
    } else {
      scores <- unlist(parallel::parLapply(cl, current_grid, cross_validation_standard,
        d = d, r = r, grid = grid, num_col = num_col, start = start, type = type,
        p = p, w = w, num_class = num_class, task = task, seed = seed, maxit_cv = maxit_cv))
      score_matrix <- matrix(scores, nrow = 1L)
    }
    scores[!is.finite(scores)] <- Inf
    list(scores = scores, score_matrix = score_matrix)
  }

  select_lambda <- function(current_grid, scores) {
    index_min <- which.min(scores); cv_min <- scores[index_min]; threshold <- cv_min + cv_tolerance
    eligible <- which(is.finite(scores) & scores <= threshold)
    selected <- eligible[which.max(current_grid[eligible])]
    list(index_min = index_min, lambda_min = current_grid[index_min], cv_min = cv_min,
         threshold = threshold, eligible = eligible, selected = selected,
         lambda_selected = current_grid[selected], cv_selected = scores[selected])
  }

  broad_lambda_grid <- as.numeric(lambda_grid)
  broad_eval <- evaluate_grid(broad_lambda_grid)
  broad_sel <- select_lambda(broad_lambda_grid, broad_eval$scores)

  if (refined_grid != 0) {
    n_grid <- length(broad_lambda_grid); j <- broad_sel$selected
    if (j == 1L) bounds <- c(broad_lambda_grid[1L] / 10, broad_lambda_grid[2L])
    else if (j == n_grid) bounds <- c(broad_lambda_grid[n_grid - 1L], broad_lambda_grid[n_grid] * 10)
    else bounds <- broad_lambda_grid[c(j - 1L, j + 1L)]
    final_grid <- exp(seq(log(bounds[1L]), log(bounds[2L]), length.out = refined_grid_length))
    final_eval <- evaluate_grid(final_grid); final_sel <- select_lambda(final_grid, final_eval$scores)
  } else {
    final_grid <- broad_lambda_grid; final_eval <- broad_eval; final_sel <- broad_sel
  }

  Estimation <- param_estim(d, r, grid, final_sel$lambda_selected, num_col, start, NULL,
                            type, p, w_total, task, seed, maxit_final)

  list(lambda_optim = final_sel$lambda_selected, selected_index = final_sel$selected,
       cv_score_selected = final_sel$cv_selected, lambda_min = final_sel$lambda_min,
       index_min = final_sel$index_min, cv_min = final_sel$cv_min,
       cv_tolerance = cv_tolerance, cv_threshold = final_sel$threshold,
       eligible_indices = final_sel$eligible, Estimation = Estimation,
       lambda_grid = final_grid, cv_scores = final_eval$scores,
       cv_score_matrix = final_eval$score_matrix,
       broad_lambda_grid = broad_lambda_grid, broad_cv_scores = broad_eval$scores,
       broad_cv_score_matrix = broad_eval$score_matrix,
       broad_index_min = broad_sel$index_min, broad_lambda_min = broad_sel$lambda_min,
       broad_cv_min = broad_sel$cv_min, broad_cv_threshold = broad_sel$threshold,
       broad_eligible_indices = broad_sel$eligible, broad_selected_index = broad_sel$selected,
       broad_lambda_selected = broad_sel$lambda_selected,
       broad_cv_score_selected = broad_sel$cv_selected,
       refined_grid_used = refined_grid != 0,
       refined_lambda_grid = final_grid, refined_cv_scores = final_eval$scores,
       refined_cv_score_matrix = final_eval$score_matrix,
       refined_index_min = final_sel$index_min, refined_lambda_min = final_sel$lambda_min,
       refined_cv_min = final_sel$cv_min, refined_selected_index = final_sel$selected)
}

# ============================================================
# 4. Increase the number of columns and return the last valid fit
# ============================================================
main_oversteps <- function(lambda_grid, N, grid, start = NULL,
                           type = c("SSR_row_HR", "SSR_row_log"), k, p,
                           num_class = 5, cl, type_noise = NULL, shape = NULL,
                           d, r, sparse_proportion, X, q, R, task, seed,
                           max_num_col = 5, zero_tolerance = 1e-8,
                           maxit_cv = 250, maxit_final = 1000,
                           cv_tolerance = 0.001, refined_grid_length = 15,
                           refined_grid = 0, verbose = TRUE) {
  type <- match.arg(type)
  w <- W_calculus(k = k, num_class = num_class, X = X, grid = grid, q = q)
  w_total <- vapply(seq_len(q), function(m) stdfEmp(R, k, grid[m, ]), numeric(1))
  diagnostics <- vector("list", max_num_col); last_valid <- NULL; last_num_col <- 0L

  for (num_col in seq_len(max_num_col)) {
    scaled_grid <- sqrt(log(d * num_col)) * lambda_grid
    fit <- main_fit(X, w, w_total, scaled_grid, grid, num_col,
                    if (num_col == 1L) start else NULL, type, k, p, num_class, cl,
                    d, r, task, seed, maxit_cv, maxit_final, cv_tolerance,
                    refined_grid_length, refined_grid)
    A_hat <- fit$Estimation$pls_matrix
    zero_columns <- colSums(abs(A_hat) > zero_tolerance) == 0L
    diagnostics[[num_col]] <- c(fit, list(num_col = num_col,
      num_parameters = d * num_col, estimated_matrix = A_hat,
      column_norms = sqrt(colSums(A_hat^2)), zero_columns = zero_columns))

    if (verbose) {
      cat("num_col =", num_col, "| lambda =", fit$lambda_optim, "\n")
      print(A_hat); cat("zero columns:", which(zero_columns), "\n")
    }
    if (any(zero_columns)) break
    last_valid <- fit; last_num_col <- num_col
  }

  diagnostics <- diagnostics[!vapply(diagnostics, is.null, logical(1))]
  if (is.null(last_valid)) { last_valid <- fit; last_num_col <- 1L }
  last_valid$num_col <- last_num_col; last_valid$diagnostics <- diagnostics
  last_valid
}
