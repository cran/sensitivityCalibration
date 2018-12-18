\name{calibrate_anim}
\alias{calibrate_anim}
\title{Make the dynamic calibration plot.}
\usage{
calibrate_anim(border, q, u, p, degree, xmax, ymax, data_matched)
}
\arguments{
\item{border}{Border or frontier of the sensitivity parameters for a fixed p.}
\item{q}{ Number of matched covariates plus treatment.}
\item{u}{ Unmeasured confounder; u = c(1,0) if the unmeasured confounder is assumed to be binary.}
\item{p}{ The probability vector corresponding to u; p = c(0.5, 0.5) if the unmeasured confounder is assumed to be Bernoulli(0.5).}
\item{degree}{Degree of freedom of the spline fit for the boundary.}
\item{xmax}{Maximum xlim of the plot.}
\item{ymax}{Maximum ylim of the plot.}
\item{data_matched}{The matched dataset.}
}

\description{
This is another main function in the package. For a given p and the border of the sensitivity parameters (lambda, delta),
a calibration plot is made for each (lambda, delta) pair on the border.
}


\details{
border is the dataframe returned by the function find_border. It has to contain at least (k+1) different
lambda/delta pairs in order to fit a smoothing spline with k dfs.
}


\examples{
\donttest{
data(NHANES_blood_lead_small_matched)
attach(NHANES_blood_lead_small_matched)

# Prepare the border
lambda_vec = c(seq(0.1,1.9,0.1), 2.2, 2.5, 3, 3.5, 4)
delta_vec = c(7.31, 5.34, 4.38, 3.76, 3.18, 2.87, 2.55, 2.36, 2.16, 1.99, 1.86,
1.74, 1.63, 1.54, 1.44, 1.40, 1.31, 1.28, 1.22, 1.08, 0.964, 0.877, 0.815, 0.750)
border = data.frame(lambda_vec, delta_vec)

calibrate_anim(border, 9, c(1,0), c(0.5,0.5), 10, 5, 3.5, NHANES_blood_lead_small_matched)

detach(NHANES_blood_lead_small_matched)
}
}
