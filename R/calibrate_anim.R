calibrate_anim <- function(border, q, u, p, degree, xmax, ymax, data_matched){

  res = NULL
  hori_line = NULL
  verti_line = NULL
  lambda_vec = border[,1]
  delta_vec = border[,2]
  for (i in 1:length(lambda_vec)) {
    # Obtain coefficients of observed coefficients
    all_coef = EM_Algorithm(q, u, p, lambda_vec[i], delta_vec[i], data_matched, all_coef = TRUE)
    coef_lambda = all_coef[[2]]
    coef_delta = all_coef[[1]]
    group = rep(i, length(coef_lambda))
    points_added = data.frame(abs(coef_lambda), abs(coef_delta), group)
    res = rbind(res, points_added)

    hori_line = rbind(hori_line, c(lambda_vec[i], delta_vec[i], 0, delta_vec[i]))
    verti_line = rbind(verti_line, c(lambda_vec[i], delta_vec[i], lambda_vec[i], 0))
  }
  res = cbind(res, rownames(res))
  colnames(res) <- c('lambda', 'delta', 'curve', 'X_obs')
  res$curve = as.factor(res$curve)

  points_added_2 <- data.frame(lambda_vec, delta_vec, seq(1,length(lambda_vec)))
  colnames(points_added_2) <- c('lambda','delta','curve')
  points_added_2$curve = as.factor(points_added_2$curve)

  hori_line <- as.data.frame(cbind(hori_line, seq(1,length(lambda_vec))))
  colnames(hori_line) <- c('x1','y1', 'x2', 'y2', 'curve')
  verti_line <- as.data.frame(cbind(verti_line, seq(1,length(lambda_vec))))
  colnames(verti_line) <- c('x1','y1', 'x2', 'y2', 'curve')

  sm = smooth.spline(border$lambda_vec, border$delta_vec, df = degree)
  border$lambda_vec = sm$x
  border$delta_vec = sm$y
  p = ggplot(border, aes(lambda_vec, delta_vec)) + geom_line(size = 1) +
    xlim(c(0, xmax)) + ylim(c(0, ymax)) +
    xlab('lambda') + ylab('delta') +
    geom_point(data = res, aes_string(x = 'lambda', y = 'delta', frame = 'curve', ids = 'X_obs'), size = 3, fontface = 'bold') +
    geom_point(data = points_added_2, aes_string(x = 'lambda', y = 'delta', frame = 'curve'), color = 'red', shape = 8, size = 3)
  p = ggplotly(p)
  return(p)
}
