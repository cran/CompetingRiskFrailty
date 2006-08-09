\name{CompetingRiskFrailtySurv}
\alias{CompetingRiskFrailtySurv}
\title{Creates a Competing Risk Model Object
}
\description{Creates an object of class 'CompetingRiskFrailtySurv' for using as a response in the model formula
in the call to the 'CompetingRiskFrailtySurvitCreate' function.
}
\usage{
CompetingRiskFrailtySurv(ID=ID,surv.time=surv.time, status=status)
}
\arguments{
  \item{ID}{identification of clusterd observations.}
  \item{surv.time}{time to follow up for right censored data.}
  \item{status}{a data frame or a matrix with columns indicating the resulting event (or comepting risk) and coded as (0,1)-indicator vectors,
                1 for event of the specified type, and 0 if such an event has not been observed.}
}
\details{The function handls only right censored data. All arguments must be 
numeric vectors. NA's can be supplied but will not be supported  by  the function 'CompetingRiskFrailtySurvfitCreate' in the further analysis.
}
\value{An object of class 'CompetingRiskFrailtySurv'.}
\references{Kauermann G. and Khomski P. (2006): \emph{Full time or part time reemployment: a competing risk model with frailties
and smooth effects using a penalty based approach}, to appear.
}
\author{Pavel Khomski <pkhomski@wiwi.uni-bielefeld.de>}
\seealso{ \code{\link{CompetingRiskFrailtySurvfitCreate}}}
\examples{
data(unemployed.comp.risk)
ID<-unemployed.comp.risk$id
surv.time<-unemployed.comp.risk$surv.time
status<-unemployed.comp.risk[,3:4]
my.Surv<-CompetingRiskFrailtySurv(ID=ID,surv.time=surv.time,status=status)
}
\keyword{misc}

