calibrate_one <- function(lambda_vec, delta_vec, q, u, p, lambda, delta, label_vec, data_matched){
  'For one particular (p, lambda, delta) triple, obtain corresponding coefficients of observed
  coefficients and plot them with the legend added.'
  border = data.frame(lambda_vec, delta_vec)

  # Obtain coefficients of observed coefficients
  all_coef = EM_Algorithm(q, u, p, lambda, delta, data_matched, all_coef = TRUE)
  coef_lambda = all_coef[[2]]
  coef_delta = all_coef[[1]]
  points_added = data.frame(abs(coef_lambda), abs(coef_delta))
  colnames(points_added) = c('lambda', 'delta')
  # print(points_added)

  # Use smooth.spline to smooth the boundary
  sm = smooth.spline(border$lambda_vec, border$delta_vec, df = 6)
  border$lambda_vec = sm$x
  border$delta_vec = sm$y

  # update delta according to the smoothing
  delta = predict(sm, x = lambda)$y

  #Plot the boundary
  border_plot = ggplot(border, aes(lambda_vec, delta_vec)) + geom_line(size = 1) +
    xlim(c(-0.5, 4.5)) + ylim(c(0, 7)) + xlab('lambda') + ylab('delta')

  # Add points to the boundary plot
  calibrate_plot = border_plot + geom_point(data = points_added, aes(x = lambda, y = delta), size = 1) +
    geom_text_repel(data = points_added, aes(lambda, delta, label = label_vec), size = 3, fontface = 'bold')

  # Add corresponding (lambda, delta) pair to the plot
  add_a_point = data.frame(lambda, delta)
  colnames(add_a_point) <- c('lambda', 'delta')
  calibrate_plot = calibrate_plot + geom_point(data = add_a_point, aes(x = lambda, y = delta), color = 'red', shape = 8, size = 3)

  # Add horiontal and vertical dotted lines
  #hori_line = data.frame('x1' = lambda, 'y1' = delta, 'x2' = -0.5, 'y2' = delta)
  #verti_line = data.frame('x1' = lambda, 'y1' = delta, 'x2' = lambda, 'y2' = 0)
  #calibrate_plot = calibrate_plot +
  #  geom_segment(aes(x = 'x1', y = 'y1', xend = 'x2', yend = 'y2'), data = hori_line, linetype = 'dotted', color = 'red', size = 1) +
  #  geom_segment(aes(x = 'x1', y = 'y1', xend = 'x2', yend = 'y2'), data = verti_line, linetype = 'dotted', color = 'red', size = 1)

  return(calibrate_plot)
}
