\name{calibrate_one}
\alias{calibrate_one}
\title{Make the calibration plot.}
\usage{
calibrate_one(lambda_vec, delta_vec, q, u, p, lambda, delta, label_vec, data_matched)
}
\arguments{
\item{lambda_vec}{A vector of lambdas that define the border.}
\item{delta_vec}{A vector of deltas that define the border.}
\item{q}{ Number of matched covariates plus treatment.}
\item{u}{ Unmeasured confounder; u = c(1,0) if the unmeasured confounder is assumed to be binary.}
\item{p}{ The probability vector corresponding to u; p = c(0.5, 0.5) if the unmeasured confounder is assumed to be Bernoulli(0.5).}
\item{lambda}{Sensitivity parameter that controls association between U and treatment assignment.}
\item{delta}{Sensitivity parameter that controls association between U and response.}
\item{label_vec}{A vector of characters of length q-1 consists of the names of observed/matched covariates.}
\item{data_matched}{The matched dataset.}
}

\description{
This is the main function in the package. Given a matched dataset and one particular (p, lambda, delta) triple,
obtain corresponding coefficients of observed coefficients and plot them with the lengend added.
This graph is meant to provide an intuitive interpretation of the magnitude of the sensitivity parameters
lambda and delta by contrasting them with the estimated coefficients of the observed covariates.
}


\details{
border is the dataframe returned by the function find_border. It has to contain at least 7 different
lambda/delta pairs in order to fit a smoothing spline with 6 dfs.

lambda and delta is a pair on the border.

label_vec is typically taken to be the columns names of the dataset, i.e., the names of the q - 1 observed covariates.
}


\examples{
data(NHANES_blood_lead_small_matched)
attach(NHANES_blood_lead_small_matched)

# Prepare the lambda_vec and delta_vec
lambda_vec = c(seq(0.1,1.9,0.1), 2.2, 2.5, 3, 3.5, 4)
delta_vec = c(7.31, 5.34, 4.38, 3.76, 3.18, 2.87, 2.55, 2.36, 2.16, 1.99, 1.86,
1.74, 1.63, 1.54, 1.44, 1.40, 1.31, 1.28, 1.22, 1.08, 0.964, 0.877, 0.815, 0.750)

calibrate_one(lambda_vec, delta_vec, 9, c(1,0), c(0.5,0.5), 1, 0.492,
colnames(NHANES_blood_lead_small_matched)[1:8], NHANES_blood_lead_small_matched)

detach(NHANES_blood_lead_small_matched)
}
