\name{CI_block_boot}
\alias{CI_block_boot}
\title{Construct the 95\% confidence interval of the treatment effect given the set of sensitivity parameters.}
\usage{
CI_block_boot(q, u, p, lambda, delta, data_matched, n_boot = 2000)
}
\arguments{
\item{q}{ Number of matched covariates plus treatment.}
\item{u}{ Unmeasured confounder; u = c(1,0) if the unmeasured confounder is assumed to be binary.}
\item{p}{ The probability vector corresponding to u; p = c(0.5, 0.5) if the unmeasured confounder
is assumed to be Bernoulli(0.5).}
\item{lambda}{Sensitivity parameter that controls association between U and treatment assignment.}
\item{delta}{Sensitivity parameter that controls association between U and response.}
\item{data_matched}{The dataset after matching.}
\item{n_boot}{Number of boostrap samples.}
}

\description{
This is the main function in the package. Given a dataset and sensitivity parameters (p, lambda, delta),
the function returns 95\% CI for the estimated treatment effect.
}


\details{
If the number of matched covariates is k, then q = k + 1.

If the hypothesized unmeasured confounder is binary, then u = c(1,0) and p = c(p, 1-p).

data_matched should be in the following format: the first (q-1) columns are matched covariates,
the qth column is the treatment status, and the (q+1)th column is the response.
See the NHANES_blood_lead_small_matched dataset for an example.

Note the input for this function is a dataset before matching. To run this function, optmatch package needs to be installed and loaded.
}


\examples{
\donttest{
data(NHANES_blood_lead_small_matched)
attach(NHANES_blood_lead_small_matched)

CI_block_boot(9, c(1,0), c(0.5,0.5), 0, 0, NHANES_blood_lead_small_matched, n_boot = 10)

detach(NHANES_blood_lead_small_matched)
}
}
