\name{NHANES_blood_lead_small_matched}
\alias{NHANES_blood_lead_small_matched}
\docType{data}
\title{
NHANES_blood_lead_small data after matching.
}
\description{
NHANES_blood_lead_small data after a full matching using the optmatch package
}
\usage{data(NHANES_blood_lead_small_matched)}
\format{
  NHANES_blood_lead_small dataset after a full matching. It consists of 500 instances and the following 12 variables:
   \describe{
   \item{\code{COP}}{treatment, 1 if cotinine level is between 0.563-14.9 ng/ml and 0 otherwise}
    \item{\code{DMARETHN}}{1 if white, 0 if others}
     \item{\code{DMPPIR}}{Poverty income ratio}
      \item{\code{HFE1}}{1 if the house is built before 1974, 0 if after 1974}
    \item{\code{HFE2}}{number of rooms in the house}
    \item{\code{HFHEDUCR}}{education level of the reference adult}
    \item{\code{HSAGEIR}}{age at the time of interview}
    \item{\code{HSFSIZER}}{size of the family}
    \item{\code{HSSEX}}{1 if male, 0 if female}
    \item{\code{PBP}}{blood lead level}
    \item{\code{U0}}{placeholder for the hypothesized unmeasured confounder U}
    \item{\code{matches}}{matched set assignment}
  }
}
\details{
We perform a full matching on the NHANES_blood_lead_small dataset using the optmatch package.
The code for constructing this matched dataset from the original dataset is given in the examples section.
We add a column U0 as placeholder for the unmeasurefor confounder U.
}

\examples{
\dontrun{
# To run this example, optmatch must be installed
set.seed(1)
library(optmatch)
data(NHANES_blood_lead_small)
attach(NHANES_blood_lead_small)

# Perform a fullmatch
fm = fullmatch(COP ~. , data = NHANES_blood_lead_small[, 1:9], min.controls = 1/4, max.controls = 4)
NHANES_blood_lead_small_matched = cbind(NHANES_blood_lead_small, matches = fm)

# Add a U0 row
U0 = rep(1, dim(NHANES_blood_lead_small_matched)[1])
NHANES_blood_lead_small_matched = cbind(NHANES_blood_lead_small_matched[,1:9], U0,
NHANES_blood_lead_small_matched[, 10:11])
}
}

\source{
NHANES III, the Third US National Health and Nutrition Examination Survey.
}
\references{
D. M. Mannino, R. Albalak, S. D. Grosse, and J. Repace. Second-hand smoke exposureand blood lead levels in U.S. children.Epidemiology, 14:719-727, 2003

A. Gelman. Scaling regression inputs by dividing by two standard deviations.Statisticsin Medicine, 27:2865-2873, 2008.
}

\keyword{datasets}
