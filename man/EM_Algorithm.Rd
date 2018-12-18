\name{EM_Algorithm}
\alias{EM_Algorithm}
\title{Estimate the treatment effect for a matched dataset given the set of sensitivity parameters.}
\usage{
EM_Algorithm(q, u, p, lambda, delta, data_matched, all_coef = FALSE,
             aug_data = FALSE, tol = 0.0001)
}
\arguments{
\item{q}{ Number of matched covariates plus treatment.}
\item{u}{ Unmeasured confounder; u = c(1,0) if the unmeasured confounder is assumed to be binary.}
\item{p}{ The probability vector corresponding to u; p = c(0.5, 0.5) if the unmeasured confounder is assumed to be Bernoulli(0.5).}
\item{lambda}{Sensitivity parameter that controls association between U and treatment assignment.}
\item{delta}{Sensitivity parameter that controls association between U and response.}
\item{data_matched}{A matched dataset. See details below.}
\item{all_coef}{TRUE then all estimated coefficients are returned, FALSE then only the estimated treatment effect is returned.}
\item{aug_data}{TRUE then the augmented dataframe at the time of convergence is returned.}
\item{tol}{Tolerance for the algorithm convergence.}
}

\description{
This is the main function in the package. Given a matched dataset and sensitivity parameters (p, lambda, delta), the function runs the EM algorithm by the method of weights and return estimated coefficients of the propensity score model and the outcome regression model.
}


\details{
If the number of matched covariates is k, then q = k + 1.

If the hypothesized unmeasured confounder is binary, then u = c(1,0) and p = c(p, 1-p).

data_matched should be in the following format: the first (q-1) columns are matched covariates, the qth column is the treatment status, the (q+1)th column is the column of unmeasured confounders U0, the (q+2)th column is the response, the last column, i.e., (q+3)th column, is the assignment of the matched set. We use the fullmatch function in the package optmatch to perform the fullmatching. See NHANES_blood_lead_small_matched for an example of a matched dataset and the examples section therein for instructions on how to construct such a matched dataset.
}


\examples{

data(NHANES_blood_lead_small_matched)
attach(NHANES_blood_lead_small_matched)

# Run the EM algorithm assuming no unmeasured confounding, i.e., lambda =delta = 0
EM_Algorithm(9, c(1,0), c(0.5,0.5), 0, 0, NHANES_blood_lead_small_matched)

# Run the EM algorithm assuming the magnitude of the unmeasured confounding is lambda =delta = 1
EM_Algorithm(9, c(1,0), c(0.5,0.5), 1, 1, NHANES_blood_lead_small_matched)

detach(NHANES_blood_lead_small_matched)
}
