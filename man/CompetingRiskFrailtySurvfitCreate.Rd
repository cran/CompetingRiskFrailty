\name{CompetingRiskFrailtySurvfitCreate}
\alias{CompetingRiskFrailtySurvfitCreate}
\title{Competing Risk Model with Frailties and Smooth Varying Coefficients
}
\description{This function offers a strategy of fitting the competing risks models with frailties for estimating a hazard function. 
The effects of covariates are modelled as varying coefficients. For flexible smoothing of the effects truncated polynomial splines are employed.
The frailties (or unobserved heterogenities) of the (clustered) observations are also modelled in a flexible way by employing a mixture of gamma
distributions. The estimation routine is based on the EM-algorithm.
}
\usage{CompetingRiskFrailtySurvfitCreate(formula=formula, data=data, na.action=na.fail,
                                         control=control, risk.names=risk.names)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details of the
                 model specification are given below.}
  \item{data}{an optional data frame containing variables in the model. If not supplied the handling of missing values will be accordingly 
              to the value of 'na.action'.}
  \item{na.action}{an optional parameter for handling with missing data.}
  \item{control}{an optional list argument with control parameters needed for optimization routine. If not specified,
                 default values from 'CompetingRiskFrailtySurvfitControl' will be taken; otherwise the values, 
		 which are defined by user, will overwrite the defaults.}
  \item{risk.names}{a character vector containing the names of competing risks. Its elements should be in the order, 
in which the columns describing the resulting states (or risks) have been supplied as a 'status' parameter
in the call of the 'CompetingRiskFrailtySurv' function.
} 		 
}
\details{
A formula argument can be specified in differentt ways. The response or left part of the formula can be either an object of class 
'CompetingRiskFrailtySurv' or it can be a call to the 'CompetingRiskFrailtySurv' function with its arguments specified in it. 
The right part may be left unspecified (without '\eqn{\sim}{~}')
or it can be of '\eqn{\sim}{~}1' type. In this case no covariates are supplied to the model, hence only baseline hazard will be adopted. 
In the case with covariates these must be supplied explicitly through their names in a standard way. In all cases the 'formula' argument 
must be supplied.

If factor variables are to be considered, their reference categories should be specified bevor applied for optimization. 

A 'control' argument can be absent, in this case the defaults from the 'CompetingRiskFrailtySurvfitControl'
will be taken; or it can be either a list with some or all components of it specified by their names with
according values, or it can be a call to the 'CompetingRiskFrailtySurvfitControl' function with parameters specified in it.

The most of computation and optimization is implemented in the "internal" function 'CompetingRiskFrailtyOptim', which is 
in the body of 'CompetingRiskFrailtySurvfitCreate'. 

The working spline bases are ones consisted from truncated polynomials. Other alternative could be B-splines, but has not
been implemented yet and might be considered in the following update of the package.

The computation procedure can take considerable time, depending on data supplied. All evaluations are written in R-code, so no external
program code like C or Fortran has been used. We advice initialy to set the number of iterations not to be too large, and eventually increase it
if needed.
}
\value{An object of class 'CompetingRiskFrailtySurvit'. It has methods for 'print' and 'plot'.
}
\references{Kauermann G. and Khomski P. (2006). \emph{Full time or part time reemployment: a competing risk model with frailties
and smooth effects using a penalty based approach}, to appear.}
\author{Pavel Khomski <pkhomski@wiwi.uni-bielefeld.de>}
\seealso{ \code{\link{CompetingRiskFrailtySurv}}, \code{\link{CompetingRiskFrailtySurvfitControl}}, \code{\link{CompetingRiskFrailtySurvfitObject}}
}
\examples{
#data(unemployed.comp.risk)
#ID<-unemployed.comp.risk$id
#surv.time<-unemployed.comp.risk$surv.time
#status<-unemployed.comp.risk[,3:4]
#my.Surv<-CompetingRiskFrailtySurv(ID=ID,surv.time=surv.time,status=status)
#my.Survfit<-CompetingRiskFrailtySurvfitCreate(my.Surv~nationality+gender+age+training+university,
#                                              data=unemployed.comp.risk,risk.names=c("full.time","part.time"))
#evaluation takes ca. 25 minuts
}
\keyword{survival}
\keyword{nonparametric}
\keyword{smooth}

