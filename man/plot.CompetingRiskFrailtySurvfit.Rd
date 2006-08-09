\name{plot.CompetingRiskFrailtySurvfit}
\alias{plot.CompetingRiskFrailtySurvfit}
\title{Plots the Smoothed Varying Coefficients}
\description{'plot.CompetingRiskFrailtySurvfit' is a proposed function for plotting of smooth components of the model,
one plot for each baseline and covariate (or factor level) effect, for each competing risk.
}
\usage{plot.CompetingRiskFrailtySurvfit(x,...)
}
\arguments{
 \item{x}{object of class 'CompetingRiskFrailtySurvfit'.}
 \item{...}{additional plot parameters.}
}
\details{All plots will be made in the device, which is specified by the user. 
One can use the object components for producing his own plots. 
}
\references{Kauermann G. and Khomski P. (2006). \emph{Full time or part time reemployment: a competing risk model with frailties
and smooth effects using a penalty based approach}, to appear.}
\author{Pavel Khomski <pkhomski@wiwi.uni-bielefeld.de>}
\seealso{\code{\link{print.CompetingRiskFrailtySurvfit}}}
\examples{
#pdf(file="myplot.pdf",h=12,w=12)
#layout(matrix(1:16,ncol=4,byrow=T),widths=c(4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4))
#par(mar=c(2,2,3,1.5))
#plot(my.Survfit)
#dev.off()
}
\keyword{hplot}

