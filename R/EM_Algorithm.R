options(warn=-1)

# q is number of covariates inclding treatment but not response
EM_Algorithm <- function(q, u, p, lambda, delta, data_matched, all_coef = FALSE, aug_data = FALSE, tol = 0.0001){

  augment_data <- function(q, u, p, delta, lambda, kappa, psi, sigma, data_matched){

    expit <- function(x) return (exp(x)/(1+exp(x)))

    # We replicate each record and attach a 1 and 0 hidden bias to them.
    k = length(p)
    data_aug = data_matched[rep(seq_len(nrow(data_matched)), each=k),]
    U = rep(u, dim(data_matched)[1])
    data_aug$U0 <- U

    # This is r_i - a_i - psi(X_i) - beta*Z - delta*U
    residual = resid(psi)
    len = dim(data_matched)[1]
    weight = numeric(len*k)
    denom = NULL

    ind_treated = which(data_aug[,q] == 1)
    X_treated = as.matrix(data_aug[ind_treated, 1:(q-1)])
    p_vec = rep(p, length(ind_treated)/k)
    u_vec = rep(u, length(ind_treated)/k)
    numerator = p_vec*expit(kappa(X_treated) + lambda*u_vec)*dnorm(residual[ind_treated]/sigma)
    for (i in seq(1, length(numerator), k)) denom = c(denom, rep(sum(numerator[i:(i+k-1)]),k))
    wt_temp_treated = numerator/denom
    weight[ind_treated] = wt_temp_treated

    ind_control = which(data_aug[,q] == 0)
    X_control = as.matrix(data_aug[ind_control, 1:(q-1)])
    p_vec = rep(p, length(ind_control)/k)
    u_vec = rep(u, length(ind_control)/k)
    denom_2 = NULL
    numerator_2 = p_vec*(1 - expit(kappa(X_control) + lambda*u_vec))*dnorm(residual[ind_control]/sigma)
    for (i in seq(1, length(numerator_2), k)) denom_2 = c(denom_2, rep(sum(numerator_2[i:(i+k-1)]),k))
    wt_temp_control = numerator_2/denom_2
    weight[ind_control] = wt_temp_control
    data_aug$weight = weight
    return(data_aug)
  }

  estimate_kappa <- function(q,lambda, data_aug) {

    ofst = lambda * data_aug$U0
    wt = data_aug$weight

    # Regress Z on all matched covariates X and treatment status Z
    model_kappa = glm(Z ~., data = data_aug[,1:q], family = 'binomial',
                      offset = ofst, weights = wt)
    coef = as.vector(coefficients(model_kappa))
    kappa_est <- function(newdata){(cbind(rep(1, nrow(newdata)),as.matrix(newdata))) %*% coef }
    return (list(kappa_est, model_kappa))
  }

  estimate_psi <- function(q, delta, data_aug){

    data_aug$Y = data_aug$Y - delta*data_aug$U0
    wt = data_aug$weight

    model_psi = lm(Y ~., data = data_aug[,c(1:q, q+2, q+3)], weights = wt)
    sigma_est = sqrt(mean((resid(model_psi)^2)))
    return (list(model_psi, sigma_est))
  }

  # Rename the treatment to Z and response to Y
  colnames(data_matched)[q] = 'Z'
  colnames(data_matched)[q+2] = 'Y'

  # Set the initial guesses
  kp_0 <- glm(Z~., data = data_matched[,1:q], family = 'binomial')
  coef = as.vector(coefficients(kp_0))
  kappa_0 <- function(newdata){(cbind(rep(1, nrow(newdata)),as.matrix(newdata))) %*% coef }

  k = length(u)
  data_aug = data_matched[rep(seq_len(nrow(data_matched)), each=k),]
  U = rep(u, dim(data_matched)[1])
  data_aug$U0 <- U

  wt = rep(p, dim(data_matched)[1])
  data_aug$Y = data_aug$Y - delta*data_aug$U0

  model_psi_0 = lm(Y ~., data = data_aug[,c(1:q,q+2,q+3)], weights = wt)
  sigma_0 = sigma(model_psi_0)

  # Set current estimates equal to the initial guesses
  current_kappa = kappa_0
  current_psi = model_psi_0
  current_sig = sigma_0
  current_beta = 0
  diff = Inf

  i = 0
  while (diff >= tol) {
    # Augment the data_matched. We compute the weights of U = 1 and U = 0
    # under current model parameters.
    data_aug = augment_data(q, u, p, delta, lambda, kappa = current_kappa, psi = current_psi, sigma = current_sig, data_matched)

    # estimate kappa
    kappa_function_model = estimate_kappa(q,lambda, data_aug)
    kappa_model = kappa_function_model[[2]]
    current_kappa = kappa_function_model[[1]]

    # estimate psi
    psi_function_model = estimate_psi(q,delta, data_aug)
    current_psi = psi_function_model[[1]]
    current_sig = psi_function_model[[2]]
    i = i + 1

    new_beta = coefficients(current_psi)[q+1]
    diff = abs(current_beta - new_beta)
    #cat('current:', current_beta, 'new:', new_beta, 'diff:', diff, '\n')
    current_beta = new_beta
  }
  if (aug_data) {
    data_aug_final = augment_data(q, u, p, delta, lambda, kappa = current_kappa, psi = current_psi, sigma = current_sig, data_matched)
    return(data_aug_final)
  }
  if (!all_coef) return(coefficients(current_psi)[q+1])
  else return(list(coefficients(current_psi)[2:q], coefficients(kappa_model)[2:q]))
}
