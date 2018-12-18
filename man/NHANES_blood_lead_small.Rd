\name{NHANES_blood_lead_small}
\alias{NHANES_blood_lead_small}
\docType{data}
\title{
A random subset of NHANES_blood_lead data.
}
\description{
A random subset of NHANES_blood_lead data for the purpose of testing.
}
\usage{data(NHANES_blood_lead_small)}
\format{
  A random sample from the NHANES_blood_lead dataset. It consists of 500 instances and the same 10 variables as the NHANES_blood_lead data.
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
  }
}
\details{
We take a 500 random sample from the NHANES_blood_lead dataset. This small dataset is primarily for the purpose of testing the algorithm.
}
\source{
NHANES III, the Third US National Health and Nutrition Examination Survey.
}
\references{
D. M. Mannino, R. Albalak, S. D. Grosse, and J. Repace. Second-hand smoke exposureand blood lead levels in U.S. children.Epidemiology, 14:719-727, 2003

A. Gelman. Scaling regression inputs by dividing by two standard deviations.Statisticsin Medicine, 27:2865-2873, 2008.
}
\examples{
data(NHANES_blood_lead_small)
}
\keyword{datasets}
