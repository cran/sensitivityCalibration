\name{find_delta}
\alias{find_delta}
\title{Estimate the maximum delta for fixed sensitivity parameters p and lambda.}
\usage{
find_delta(q, u, p, lambda, start_value_low, start_value_high,
data_matched, n_boot = 200, tol = 0.01)
}
\arguments{
\item{q}{ Number of matched covariates plus treatment.}
\item{u}{ Unmeasured confounder; u = c(1,0) if the unmeasured confounder is assumed to be binary.}
\item{p}{ The probability vector corresponding to u; p = c(0.5, 0.5) if the unmeasured confounder
is assumed to be Bernoulli(0.5).}
\item{lambda}{A lambda value.}
\item{start_value_low}{Starting value for the binary search (the lower endpoint).}
\item{start_value_high}{Starting value for the binary search (the higher endpoint).}
\item{data_matched}{The dataset after matching.}
\item{n_boot}{Number of boostrap samples used to approximate the CI.}
\item{tol}{Tolerance for the binary search.}
}

\description{
Estimate the maximum delta value for a given p and lambda, so that the estimated 95\% confidence
interval for the treatment effect is still significant. Note in order to run this function,
optmatch package needs to be installed and loaded.
}


\details{
start_value_low and start_value_high are user supplied numbers to start the binary search.
}


\examples{
\donttest{
data(NHANES_blood_lead_small_matched)
attach(NHANES_blood_lead_small_matched)

find_delta(9, c(1,0), c(0.5,0.5), 1, 1, 3,
NHANES_blood_lead_small_matched, n_boot = 1000)

detach(NHANES_blood_lead_small_matched)
}
}
