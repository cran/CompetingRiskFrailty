\name{print.CompetingRiskFrailtySurvfit}
\alias{print.CompetingRiskFrailtySurvfit}
\title{Prints an Object of class 'CompetingRiskFrailtySurvfit'}
\description{Prints estimates of fixed parameters of smoothing varying coefficients, penalty values for their random parts 
as well as mixture weights of frailty mixture distribution, their degrees of feedom and optimal AIC- and marginal log likelihood value of the model.
}
\usage{print.CompetingRiskFrailtySurvfit(x,...)
}
\arguments{
  \item{x}{  object of class 'CompetingRiskFrailtySurvfit'.}
  \item{...}{  additional parameters for print.}
}
\references{Kauermann G. and Khomski P. (2006). \emph{Full time or part time reemployment: a competing risk model with frailties
and smooth effects using a penalty based approach}, to appear.}
\author{Pavel Khomski <pkhomski@wiwi.uni-bielefeld.de>}
\seealso{\code{\link{plot.CompetingRiskFrailtySurvfit}}}
\examples{
#print.CompetingRiskFrailtySurvfit(my.Survfit)
}
\keyword{print}

