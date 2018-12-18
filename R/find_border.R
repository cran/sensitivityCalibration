find_border <- function(q, u, p, lambda_vec, start_value_low, start_value_high, data_matched, n_boot = 2000, tol = 0.01){
  'For a vector of lambda values, find_border computes the maximum delta as described by the function
  find_delta. The function returns a rough border.'
  delta_vec = numeric(length(lambda_vec))
  for (i in 1:length(lambda_vec)){
    lambda = lambda_vec[i]
    delta_vec[i] = find_delta(q, u, p, lambda, start_value_low, start_value_high, data_matched, n_boot, tol)
    start_value_high = delta_vec[i]
    start_value_low = 0
  }
  return(data.frame(lambda_vec, delta_vec))
}
