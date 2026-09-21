source("C:/Github/Penalized_least-squares_estimator/R_functions/helper_functions.R")
#####Evaluate EDS when sampling and fitting mixture logistic model
branch_path <- "C:/Results_Github/Penalized_least_squares_estimator/Extreme_directions_identificcation/sample_fit_mix_log"

n <- 100
results <- vector("list", n)
score_EDS <- 0
for (i in seq_len(n)) {
  file_path <- file.path(branch_path, paste0("result_", i, ".rds"))
  res <- readRDS(file_path)
  true_matrix <- res$True_matrixA
  estim_matrix <- res$Estimation[[1]]$Estimation$pls_matrix
  true_extreme_directions <- extract_signatures(true_matrix)
  estim_extreme_directions <- extract_signatures(estim_matrix)
  EDSi <- EDS(true_extreme_directions , estim_extreme_directions)
  score_EDS <- score_EDS + EDSi/n
}
