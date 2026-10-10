branch_path <- "C:/Results_Github/Penalized_least_squares_estimator/Extreme_directions_identificcation/DA_mix_HR/extreme_direction_identification"

n <- 100
N <- seq(500, 5000, by = 500)
score_EDS_HR <- numeric(length(N))

for (i in seq_len(n)) {
  file_path <- file.path(branch_path, paste0("result_", i, ".rds"))
  res <- readRDS(file_path)

  true_matrix <- res$True_matrixA
  true_extreme_directions <- extract_signatures(true_matrix)

  for (j in seq_along(N)) {
    N_name <- paste0("N_", N[j])

    # First and only kn result
    fit <- res$Estimation[[N_name]]$Estimation[[1]]

    estim_matrix <- fit$Estimation$pls_matrix
    lambda <- fit$lambda_optim
    lambda_grid <- fit$lambda_grid

    estim_extreme_directions <- extract_signatures(estim_matrix)

    EDSij <- EDS(
      true_extreme_directions,
      estim_extreme_directions
    )

    cat(
      "For seed", i,
      "and N =", N[j],
      "the EDS is", EDSij, "\n"
    )

    score_EDS_HR[j] <- score_EDS[j] + EDSij / n
  }
}


####Evaluation of ED-S for mix-logistic

branch_path <- "C:/Results_Github/Penalized_least_squares_estimator/Extreme_directions_identificcation/DA_mix_log/Frechet_Noise_Night"

n <- 100
N <- seq(500, 5000, by = 500)
score_EDS_log <- numeric(length(N))

for (i in seq_len(n)) {
  file_path <- file.path(branch_path, paste0("result_", i, ".rds"))
  res <- readRDS(file_path)

  true_matrix <- res$True_matrixA
  true_extreme_directions <- extract_signatures(true_matrix)

  for (j in seq_along(N)) {
    N_name <- paste0("N_", N[j])

    # First and only kn result
    fit <- res$Estimation[[N_name]]$Estimation[[1]]

    estim_matrix <- fit$Estimation$pls_matrix
    lambda <- fit$lambda_optim
    lambda_grid <- fit$lambda_grid

    estim_extreme_directions <- extract_signatures(estim_matrix)

    EDSij <- EDS(
      true_extreme_directions,
      estim_extreme_directions
    )

    cat(
      "For seed", i,
      "and N =", N[j],
      "the EDS is", EDSij, "\n"
    )

    score_EDS_log[j] <- score_EDS_log[j] + EDSij / n
  }
}

##########Plot
colors <- c("#56B4E9", "#E69F00")
x_values <- seq(500, 5000, by = 500); x_pos <- seq_along(x_values)
par(mgp = c(2, 0.7, 0))
plot(x_pos, score_EDS_log, type = "o", col = colors[1], pch = 16, lwd = 2.5, cex = 1.5, ylim = c(0, 1), xlab = "Sample size n", ylab = "ED-S", cex.lab = 1.5, cex.axis = 1.3, xaxt = "n")
axis(1, at = which(x_values %% 1000 == 0), labels = x_values[x_values %% 1000 == 0], cex.axis = 1.3)
lines(x_pos, score_EDS_HR, type = "o", col = colors[2], pch = 16, lwd = 2.5, cex = 1.5)

