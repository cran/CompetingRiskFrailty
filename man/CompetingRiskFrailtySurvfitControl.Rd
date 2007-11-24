\name{CompetingRiskFrailtySurvfitControl}
\alias{CompetingRiskFrailtySurvfitControl}
\title{Control Values for Fitting of Competing Risks with Frailties Model}
\description{The values supplied in the call of 'CompetingRiskFrailtySurvfitCreate' function will replace the defaults, and
a list with all possible arguments is returned. The returned list
is used as the 'control' argument to the 'CompetingRiskFrailtySurvfitCreate' function.
}
\usage{
CompetingRiskFrailtySurvfitControl(niter.EM=50, niter.epoch = 2, tol.epoch = 1e-08, tol.variance = 1e-08,
                                   tol.frailty = 1e-06, print.penalty.mixture=TRUE,print.EM=TRUE,
				   print.estimates=FALSE,print.log.lik=TRUE,...)
}
\arguments{
  \item{niter.EM}{maximum number of the (outer) EM-iterations.}
  \item{niter.epoch}{maximum number of the (inner) iterations in optimization for 
                     varying coefficients parameters theta and penalty parameters of their random parts, within an EM-iteration.}		     
  \item{tol.epoch}{tolerance for the convergence criterion for the fixed and random parameters of the varying coefficients.}
  \item{tol.variance}{tolerance for the convergence criterion for the penalty values of varying coefficients.}
  \item{tol.frailty}{tolerance for the the convergence criterion for the frailty terms.}
  \item{print.penalty.mixture}{logical value for printing the value of the penalty parameter from the specified grid of values.}
  \item{print.EM}{logical value for printing the current number of an EM-itertion.}
  \item{print.estimates}{logical value for printing the estimates of the fixed parameters theta and penalties
                         of varying coefficients after the last EM-iteration.}
  \item{print.log.lik}{logical value for printing the marginal log likelihood of the model in each EM-iteration.}
  \item{\dots}{other parameters which can only be: 'num.knots' for the number of spline knots for 
               survival time. If specified it is a vector of integer values of the length equal to the number of competing risks.
	       If not specified, the optimal values will be defined internally.
  }
}
\details{The defaults or user specified values are applied as the 'control' argument in the
call of the 'CompetingRiskFrailtySurvfitCreate' function. It can be an (empty) list object or a call to the 
'CompetingRiskFrailtySurvfitControl' function itself, whether or not with supplied arguments to be changed from their default values.
}
\value{a list with components for each of the possible arguments.
}
\references{Kauerman G. and Khomski P. (2006). \emph{Full time or part time reemployment: a competing risk model with frailties
and smooth effects using a penalty based approach}, to appear.}
\author{Pavel Khomski <pkhomski@wiwi.uni-bielefeld.de>}
\seealso{\code{\link{CompetingRiskFrailtySurvfitCreate}}}
\keyword{list}

