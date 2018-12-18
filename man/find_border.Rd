\name{find_border}
\alias{find_border}
\title{Find the lambda-delta boundary for a fixed sensitivity parameter p.}
\usage{
find_border(q, u, p, lambda_vec, start_value_low, start_value_high,
data_matched, n_boot = 2000, tol = 0.01)
}
\arguments{
\item{q}{ Number of matched covariates plus treatment.}
\item{u}{ Unmeasured confounder; u = c(1,0) if the unmeasured confounder is assumed to be binary.}
\item{p}{ The probability vector corresponding to u; p = c(0.5, 0.5) if the unmeasured confounder
is assumed to be Bernoulli(0.5).}
\item{lambda_vec}{A sequence of lambda values.}
\item{start_value_low}{Starting value for the binary search (the lower endpoint).}
\item{start_value_high}{Starting value for the binary search (the higher endpoint).}
\item{data_matched}{The dataset after matching.}
\item{n_boot}{Number of boostrap samples used to approximate the CI.}
\item{tol}{Tolerance for the binary search.}
}

\description{
Given the dataset, unmeasured confounder, sensitivity parameter p, and a sequence of lambda values,
the function uses binary search to find a sequence of delta corresponding to each lambda
in the lambda_vec such that the estimated 95\% for the treatment effect barely covers 0.
The function returns a dataframe consisting of lambda_vec and the corresponding deltas.
See below for an example.
}


\details{
start_value_low and start_value_high are user supplied numbers to start the binary search.
}

\examples{
\donttest{
data(NHANES_blood_lead_small_matched)
attach(NHANES_blood_lead_small_matched)

find_border(9, c(1,0), c(0.5,0.5), c(0.5,1,1.5), 0, 4,
NHANES_blood_lead_small_matched, n_boot = 1000)

detach(NHANES_blood_lead_small_matched)
}
}
