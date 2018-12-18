CI_block_boot <- function(q, u, p, lambda, delta, data_matched, n_boot = 2000){
  # match the dataset
  beta_boot = numeric(n_boot)

  for (i in 1:n_boot){
    #cat('i',i,'\n')
    # resample the matched sets
    ind = sample(unique(data_matched$matches), length(unique(data_matched$matches)), replace = TRUE)
    new_data = data.frame()
    for (id in ind){
      # If the matched set has been sampled, generate a random id for this new one
      if (id %in% new_data$matches) {
        sub_data_frame = data_matched[data_matched$matches == id,]
        sub_data_frame$matches = stri_rand_strings(1,4)
        new_data = rbind(new_data, sub_data_frame)
      } else
        new_data = rbind(new_data, data_matched[data_matched$matches == id,])
    }

    beta_boot[i] = EM_Algorithm(q, u, p, lambda, delta, new_data)
  }
  CI = quantile(beta_boot, c(0.025, 0.975))
  return(CI)
}
