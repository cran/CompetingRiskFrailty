\name{CompetingRiskFrailtySurvfitObject}
\alias{CompetingRiskFrailtySurvfitObject}
\title{Fitted Object}
\description{An object of the class 'CompetingRiskFrailtySurvfit' returned as a list by the 'CompetingRiskFrailtySurvfitCreate' function
and representing a fitted competing risks hazard model with frailties and varying coefficients. 
It has methods for generic functions print and plot.
}
\value{

  \item{L}{number of comepting risks.}
  \item{M.list}{a list with elements specifying the number of chosen gamma distributions of frailties for each competing risk.}
  \item{fixed.coef.optim}{named vector of estimates of fixed parameters, beta, of the model.}
  \item{random.coef.optim}{named vector of estimates of random parameters, u, of the model.}
  \item{penalty.varying.optim}{named vector of penalty values for the random parts of the varying coefficients.}
  \item{penalty.weights.optim}{optimal value of the penalty parameter of the frailty mixture density.}
  \item{grid.frame}{a frame of grid values for plotting of the varying coefficients.}
  \item{varying.list}{a frame containing the values of varying coefficients.}
  \item{deviation.list}{a list containing the values of deviation for varying coefficients.}
  \item{frailty.list}{a list of frailty estimates, each component is a vector according to the competing risk.}
  \item{mixture.weights}{a vector of estimated mixture weights for the components of the mixture frailty distribution.}
  \item{aic.optim}{an optimal aic value of the model.}
  \item{df.weights.optim}{an optimal value of the degrees of freedom for mixture weights.}
  \item{log.lik.margin.optim}{an optimal value of the marginal log likelihood of the model.}
  \item{p}{number of covariates; for factors including their categories (excluding reference category).}
  \item{factor.names}{covariates' names; for factors the names of categories (excluding reference category).}
  \item{risk.names}{a character vector with the names of competing risks.}
  
  } 
\references{Kauermann G. and Khomski P. (2006). \emph{Full time or part time reemployment: a competing risk model with frailties
and smooth effects using a penalty based approach} , to appear.
}
\author{Pavel Khomski <pkhomski@wiwi.uni-bielefeld.de>}
\seealso{\code{\link{CompetingRiskFrailtySurv}}, \code{\link{CompetingRiskFrailtySurvfitCreate}}}
\keyword{misc}  

