\name{NHANES_blood_lead}
\alias{NHANES_blood_lead}
\docType{data}
\title{
Second hand smoking and blood lead levels dataset from NHANES III.
}
\description{
A dataset constructed from NHANES III.
}
\usage{data(NHANES_blood_lead)}
\format{
  A data frame with 4519 observations on the following 10 variables.
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
We follow Mannino rt al. (2003) in constructing a dataset that includes children aged 4-16 years old for whom both serum cotinine levels and blood lead levels were measured in the Third National Health and Nutrition Examination Survey (NHANES III), along with the following variables: race/ethnicity, age, sex, poverty income ratio, education level of the reference adult, family size, number of rooms in the house, and year the house was constructed. The biomarker cotinine is a metabolite of nicotine and an indicator of second-hand smoke exposure. Treatment status is 1 if cotinine level is between 0.563-14.9 ng/ml and 0 otherwise. All continuous/ordinal variables are standardized by subtracting the mean and divided by 2 standard deviations so that they are more comparable to binary covariates (Gelman 2008).
}
\source{
NHANES III, the Third US National Health and Nutrition Examination Survey.
}
\references{
D. M. Mannino, R. Albalak, S. D. Grosse, and J. Repace. Second-hand smoke exposureand blood lead levels in U.S. children.Epidemiology, 14:719-727, 2003

A. Gelman. Scaling regression inputs by dividing by two standard deviations.Statisticsin Medicine, 27:2865-2873, 2008.
}
\examples{
data(NHANES_blood_lead)
}
\keyword{datasets}
