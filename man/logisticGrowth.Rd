% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sourceSHAPE.R
\name{logisticGrowth}
\alias{logisticGrowth}
\title{This function is simply an implementation of the logistic growth equation where:
f(x) = K / (1 + ((K - N_0)/N_0) *exp-k(x-x_0))  ; Where x_0 is an adjustment to the position of the midpoint of the curve's maximum value
K = the curves maximum value, k = the steepness of the curve (growth rate), and N_0 is the starting population
it includes parameters to change the midpoint as well as change the natural exponent (ie: exp) to some other value.
NOTE: This is for continuous growth, and since SHAPE is discrete at present this is an unused function.}
\usage{
logisticGrowth(func_rate, func_step, func_startPop = NULL,
  func_maxPop = NULL, func_midAdjust = 0,
  func_basalExponent = exp(1))
}
\arguments{
\item{func_rate}{The basal growth rate of individuals in the SHAPE run.}

\item{func_step}{This is the number of steps forward for which you wish to calculate the growth expected.}

\item{func_startPop}{The sum of the populations in the evolving community.}

\item{func_maxPop}{The carrying capacity of the enviromment being simulated.}

\item{func_midAdjust}{The midpoint which controls the point of inflection for the logistic equation.  Beware, change this at your own risk
as its impact will varrying based on the population sizes being simulated.  Ideally, don't change this value from its default.}

\item{func_basalExponent}{This defaults as the natural exponent "e" / "exp".  Change it at your own risk.}
}
\value{
Returns a single value representing the amount of logistic growth expected by the community
}
\description{
This function is simply an implementation of the logistic growth equation where:
f(x) = K / (1 + ((K - N_0)/N_0) *exp-k(x-x_0))  ; Where x_0 is an adjustment to the position of the midpoint of the curve's maximum value
K = the curves maximum value, k = the steepness of the curve (growth rate), and N_0 is the starting population
it includes parameters to change the midpoint as well as change the natural exponent (ie: exp) to some other value.
NOTE: This is for continuous growth, and since SHAPE is discrete at present this is an unused function.
}
\examples{
# This calculates logistic growth based on the mathematical continuous time algorithm
logisticGrowth(func_rate = 2, func_step = 1, func_startPop = 1e2, func_maxPop = 1e4)
# It normally takes log2(D) steps for a binary fission population to reach carrying capacity,
# where D is max/start, in this case D = 100 and so it should take ~ 6.64 turns
logisticGrowth(func_rate = 2, func_step = c(1,2,3,6,6.64,7), func_startPop = 1e2, func_maxPop = 1e4)

}
