# ============================================================
# Estimation, cross-validation, and overstep selection
# Required beforehand: source("helper_functions.R"), load main.dll,
# and define/load stdfEmp().
# ============================================================

param_estim_path_fold2 <- function(d, r, grid, lambda, num_col = NULL, start, start_dep = NULL,
                                   type = c("SSR_row_HR", "SSR_row_log"), p, w, task = NULL,
                                   seed = NULL, maxit = 1000) {
  type <- match.arg(type)
  if (is.null(num_col)) num_col <- r
  if (!is.null(seed)) set.seed(seed)

  q <- nrow(grid)
  l <- d * num_col
  v <- if (type == "SSR_row_HR") d * (d - 1L) / 2L else 1L

  if (length(start) != l) stop("start must have length d * num_col.")
  if (length(w) != q) stop("w must have length nrow(grid).")

  if (is.null(start_dep)) start_dep <- if (type == "SSR_row_HR") rep(0.2, v) else runif(1, 0.15, 0.25)
  if (length(start_dep) == 1L && v > 1L) start_dep <- rep(start_dep, v)
  if (length(start_dep) != v) stop("start_dep has the wrong length.")

  start <- pmin(pmax(start, 0), 1)
  start_dep <- if (type == "SSR_row_HR") pmin(pmax(start_dep, 0.1), 2) else pmin(pmax(start_dep, 0.1), 0.3)

  p_C <- as.double(p)
  lambda_C <- as.double(lambda)
  d_C <- as.integer(d)
  num_col_C <- as.integer(num_col)
  q_C <- as.integer(q)
  w_C <- as.double(w)
  grid_C <- as.double(t(grid))

  call_ssr <- function(theta_A, theta_dep, penalty = lambda_C) {
    dep_C <- if (type == "SSR_row_log") rep(theta_dep, num_col) else theta_dep

    value <- .C(type, p_C, as.double(penalty), as.double(theta_A), d_C,
                num_col_C, q_C, as.double(dep_C), w_C, grid_C,
                R = double(1))$R

    if (is.finite(value)) value else 1e16
  }

  lower_dep <- rep(0.1, v)
  upper_dep <- if (type == "SSR_row_log") rep(0.3, v) else rep(2, v)

  # Stage 1: jointly estimate A and the dependence parameters
  objective_joint <- function(theta) {
    theta_A <- theta[seq_len(l)]
    theta_dep <- theta[l + seq_len(v)]
    call_ssr(theta_A, theta_dep)
  }

  fit_joint <- optim(c(start, start_dep), objective_joint, method = "L-BFGS-B",
                     lower = c(rep(0, l), lower_dep),
                     upper = c(rep(1, l), upper_dep),
                     control = list(maxit = maxit))

  theta_A <- fit_joint$par[seq_len(l)]
  theta_dep <- fit_joint$par[l + seq_len(v)]
  A_joint <- matrix(theta_A, nrow = d, ncol = num_col, byrow = TRUE)
  dep_joint <- if (type == "SSR_row_HR") construct_symmetric_matrix_2(theta_dep) else theta_dep

  # Stop after joint estimation for extreme-direction identification
  if (identical(task, "ED_identification")) {
    return(list(pls_matrix = A_joint, pls_dep = dep_joint,
                pls_dep_vector = theta_dep, par = fit_joint$par,
                convergence = fit_joint$convergence,
                objective = fit_joint$value, counts = fit_joint$counts,
                message = fit_joint$message))
  }

  # Stage 2: fix dependence parameters and estimate A
  objective_A <- function(theta) call_ssr(theta, theta_dep)

  fit_A <- optim(theta_A, objective_A, method = "L-BFGS-B",
                 lower = rep(0, l), upper = rep(1, l),
                 control = list(maxit = maxit))

  A_final <- matrix(normalize_group(fit_A$par, num_col),
                    nrow = d, ncol = num_col, byrow = TRUE)

  A_vector <- as.double(t(A_final))

  # Stage 3: fix A and estimate dependence parameters, with penalty equal to 0
  objective_dep <- function(theta) call_ssr(A_vector, theta, penalty = 0)

  fit_dep <- optim(theta_dep, objective_dep, method = "L-BFGS-B",
                   lower = lower_dep, upper = upper_dep,
                   control = list(maxit = maxit))

  dep_vector <- fit_dep$par
  dep_final <- if (type == "SSR_row_HR") construct_symmetric_matrix_2(dep_vector) else dep_vector

  list(pls_matrix = A_final, pls_dep = dep_final,
       pls_dep_vector = dep_vector, par = c(A_vector, dep_vector),
       convergence = fit_dep$convergence, objective = fit_dep$value,
       counts = fit_dep$counts, message = fit_dep$message,
       joint_fit = fit_joint, A_fit = fit_A, dependence_fit = fit_dep)
}



cross_validation_standard <- function(lambda, d, r, grid, num_col = NULL, start,
                                      type = c("SSR_row_HR", "SSR_row_log"), p, w,
                                      num_class = 5, task, seed, maxit_cv = 250) {
  type <- match.arg(type)
  if (is.null(num_col)) num_col <- r
  q <- nrow(grid); scores <- rep(Inf, num_class)
  for (class_k in seq_len(num_class)) {
    fit <- tryCatch(param_estim(d, r, grid, lambda, num_col, start, NULL, type, p,
                                w$train[[class_k]], task, seed, maxit_cv), error = function(e) NULL)
    if (is.null(fit)) next
    dep_C <- if (type == "SSR_row_log") rep(fit$pls_dep_vector, num_col) else fit$pls_dep_vector
    scores[class_k] <- .C(type, as.double(p), as.double(0), as.double(t(fit$pls_matrix)),
                          as.integer(d), as.integer(num_col), as.integer(q), as.double(dep_C),
                          as.double(w$test[[class_k]]), as.double(t(grid)), R = double(1))$R
  }
  if (all(!is.finite(scores))) Inf else mean(scores[is.finite(scores)])
}



cross_validation_path_fold2 <- function(class_k, lambda_grid, d, r, grid, num_col = NULL, start, type = c("SSR_row_HR", "SSR_row_log"), p, w_train, w_test, task, seed, maxit_cv = 250) {
  type <- match.arg(type)
  if (is.null(num_col)) num_col <- r

  q <- nrow(grid)
  l <- d * num_col
  v <- if (type == "SSR_row_HR") d * (d - 1L) / 2L else 1L

  lambda_order <- order(lambda_grid, decreasing = TRUE)
  lambda_sorted <- lambda_grid[lambda_order]
  scores_sorted <- rep(Inf, length(lambda_sorted))
  current_A <- start
  current_dep <- if (type == "SSR_row_HR") rep(0.2, v) else 0.2

  p_C <- as.double(p)
  d_C <- as.integer(d)
  num_col_C <- as.integer(num_col)
  q_C <- as.integer(q)
  grid_C <- as.double(t(grid))
  test_w_C <- as.double(w_test[[class_k]])

  for (j in seq_along(lambda_sorted)) {
    fit <- tryCatch(
      param_estim_path_fold2(d, r, grid, lambda_sorted[j], num_col, current_A, current_dep, type, p, w_train[[class_k]], task, seed, maxit_cv),
      error = function(e) {
        message("Fold ", class_k, ", lambda ", lambda_sorted[j], ": ", conditionMessage(e))
        NULL
      }
    )

    if (is.null(fit)) next

    dep_C <- if (type == "SSR_row_log") rep(fit$pls_dep_vector, num_col) else fit$pls_dep_vector

    score <- .C(type, p_C, as.double(0), as.double(t(fit$pls_matrix)), d_C, num_col_C, q_C, as.double(dep_C), test_w_C, grid_C, R = double(1))$R

    if (is.finite(score)) scores_sorted[j] <- score

    if (length(fit$par) == l + v && all(is.finite(fit$par))) {
      current_A <- fit$par[seq_len(l)]
      current_dep <- fit$par[l + seq_len(v)]
    }
  }

  scores <- rep(Inf, length(lambda_grid))
  scores[lambda_order] <- scores_sorted
  scores
}

main_fit <- function(X, w, w_total, lambda_grid, grid, num_col = NULL, start = NULL,
                     type = c("SSR_row_HR", "SSR_row_log"), k, p, num_class = 5, cl,
                     d, r, task, seed, maxit_cv = 250, maxit_final = 1000,
                     cv_tolerance = 0.001, refined_grid_length = 15, refined_grid = 0 , type_CV = c("cross_validation_path_fold2" , "cross_validation_standard")) {
  type <- match.arg(type)
  if (is.null(num_col)) num_col <- r
  if (is.null(start)) start <- as.double(t(starting_point(X, num_col)))
  if (length(start) != d * num_col) stop("start must have length d * num_col.")
  type_CV <- match.arg(type_CV)
  evaluate_grid <- function(current_grid) {
    if (type_CV == "cross_validation_path_fold2") {
      fold_scores <- parallel::parLapply(
        cl, seq_len(num_class), cross_validation_path_fold2,
        lambda_grid = current_grid, d = d, r = r, grid = grid,
        num_col = num_col, start = start, type = type, p = p,
        w_train = w$train, w_test = w$test, task = task,
        seed = seed, maxit_cv = maxit_cv
      )

      score_matrix <- do.call(rbind, fold_scores)
      scores <- colMeans(score_matrix)

    } else {
      scores <- unlist(parallel::parLapply(
        cl, current_grid, cross_validation_standard,
        d = d, r = r, grid = grid, num_col = num_col,
        start = start, type = type, p = p, w = w,
        num_class = num_class, task = task, seed = seed,
        maxit_cv = maxit_cv
      ))

      score_matrix <- matrix(scores, nrow = 1L)
    }

    scores[!is.finite(scores)] <- Inf
    list(scores = scores, score_matrix = score_matrix)
  }


  select_lambda <- function(current_grid, scores) {
    finite <- which(is.finite(scores)); if (!length(finite)) stop("All CV scores are non-finite.")
    index_min <- finite[which.min(scores[finite])]; cv_min <- scores[index_min]
    threshold <- cv_min + cv_tolerance; eligible <- which(is.finite(scores) & scores <= threshold)
    selected <- eligible[which.max(current_grid[eligible])]
    list(index_min = index_min, lambda_min = current_grid[index_min], cv_min = cv_min,
         threshold = threshold, eligible = eligible, selected = selected,
         lambda_selected = current_grid[selected], cv_selected = scores[selected])
  }

  broad_grid <- as.numeric(lambda_grid); broad_eval <- evaluate_grid(broad_grid)
  broad_sel <- select_lambda(broad_grid, broad_eval$scores)
  if (refined_grid != 0) {
    j <- broad_sel$selected; n_grid <- length(broad_grid)
    bounds <- if (j == 1L) c(broad_grid[1L] / 10, broad_grid[2L]) else if (j == n_grid) c(broad_grid[n_grid - 1L], broad_grid[n_grid] * 10) else broad_grid[c(j - 1L, j + 1L)]
    final_grid <- exp(seq(log(min(bounds)), log(max(bounds)), length.out = refined_grid_length))
    final_eval <- evaluate_grid(final_grid); final_sel <- select_lambda(final_grid, final_eval$scores)
  } else {
    final_grid <- broad_grid; final_eval <- broad_eval; final_sel <- broad_sel
  }

  estimation <- param_estim_path_fold2(d, r, grid, final_sel$lambda_selected, num_col, start, NULL,
                                       type, p, w_total, task, seed, maxit_final)
  list(lambda_optim = final_sel$lambda_selected, selected_index = final_sel$selected,
       cv_score_selected = final_sel$cv_selected, lambda_min = final_sel$lambda_min,
       index_min = final_sel$index_min, cv_min = final_sel$cv_min,
       cv_tolerance = cv_tolerance, cv_threshold = final_sel$threshold,
       eligible_indices = final_sel$eligible, Estimation = estimation,
       lambda_grid = final_grid, cv_scores = final_eval$scores,
       cv_score_matrix = final_eval$score_matrix, broad_lambda_grid = broad_grid,
       broad_cv_scores = broad_eval$scores, broad_cv_score_matrix = broad_eval$score_matrix,
       broad_selection = broad_sel, refined_grid_used = refined_grid != 0)
}

main_oversteps <- function(lambda_grid, N, grid, start = NULL,
                           type = c("SSR_row_HR", "SSR_row_log"), k, p,
                           num_class = 5, cl, type_noise = NULL, shape = NULL,
                           d, r, sparse_proportion, X, q = nrow(grid), R, task, seed,
                           max_num_col = 5, zero_tolerance = 1e-8,
                           maxit_cv = 250, maxit_final = 1000,
                           cv_tolerance = 0.001, refined_grid_length = 15,
                           refined_grid = 0, verbose = TRUE) {
  type <- match.arg(type)
  if (q != nrow(grid)) stop("q must equal nrow(grid).")
  w <- W_calculus(k, num_class, X, grid, q)
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
      num_parameters = d * num_col + if (type == "SSR_row_HR") d * (d - 1L) / 2L else 1L,
      estimated_matrix = A_hat, column_norms = sqrt(colSums(A_hat^2)), zero_columns = zero_columns))
    if (verbose) {
      cat("num_col =", num_col, "| lambda =", fit$lambda_optim, "\n"); print(A_hat)
      cat("zero columns:", paste(which(zero_columns), collapse = ", "), "\n")
    }
    if (any(zero_columns)) break
    last_valid <- fit; last_num_col <- num_col
  }

  diagnostics <- diagnostics[!vapply(diagnostics, is.null, logical(1))]
  if (is.null(last_valid)) { last_valid <- fit; last_num_col <- 1L }
  last_valid$num_col <- last_num_col; last_valid$diagnostics <- diagnostics; last_valid
}

setup_estimation_workers <- function(cl, main_file, helper_file, dll_file) {
  paths <- normalizePath(c(main_file, helper_file, dll_file), mustWork = TRUE)
  parallel::clusterExport(cl, "paths", envir = environment())
  parallel::clusterEvalQ(cl, {
    source(paths[2L]); source(paths[1L]); library(mvtnorm); dyn.load(paths[3L]); NULL
  })
  invisible(parallel::clusterEvalQ(cl, c(HR = is.loaded("SSR_row_HR"), log = is.loaded("SSR_row_log"))))
}
