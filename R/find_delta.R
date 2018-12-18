find_delta <- function(q, u, p, lambda, start_value_low, start_value_high, data_matched, n_boot = 200, tol = 0.01){
  'For the user specified lambda, find the largest delta such that the treatment effect
  is still significant via binary search.'
  L_delta = start_value_low
  R_delta = start_value_high
  delta_final = L_delta
  while (R_delta - L_delta > tol) {
    mid_point_delta = (L_delta + R_delta)/2
    CI = CI_block_boot(q, u, p, lambda, mid_point_delta, data_matched, n_boot = n_boot)
    cat('current_lambda, current_delta', c(lambda, L_delta, R_delta, delta_final, CI), '\n')
    CI_low = quantile(CI, 0.025)
    CI_high = quantile(CI, 0.975)
    if (CI_low >= 0) {
      L_delta = mid_point_delta
      delta_final = L_delta
    } else{
      R_delta = mid_point_delta
    }
  }
  return(delta_final)
}
