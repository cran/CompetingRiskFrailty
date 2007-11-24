\name{unemployed.comp.risk}
\docType{data}
\alias{unemployed.comp.risk}
\title{Sample of German Unemployed}
\description{This data give the time spent in the unemployment spell and two competing risks: full time and part time reemployment
as possible resulting states for leaving that spell for a sample of german unemployed individuals. The full data have been recieved from the 
\emph{Deutsches Institut fuer Wirtschaftsforschung}, s. \url{www.diw.de} and were analysed in the larger study. 
}
\usage{data(unemployed.comp.risk)
}
\format{A data frame containing 500 sampled observations. Description of columns:
\itemize{
\item{id}{Observations or spells are clustered according to the identification number of the person.}
\item{surv.time}{Time spent in the unemployment spell, measured in months. The spells which lasted more than 36 months have been
 truncated on 36 months and denoted as sensored.}
\item{full.time}{binary indicator for leaving an unemloyment spell rsulting in a full time job.}
\item{part.time}{binary indicator for leaving an unemloyment spell rsulting in a part time job.}
\item{nationality}{Nationality of the unemployed person, german vs. others.}
\item{gender}{Gender of the person.}
\item{age}{Age has been categorized in young (till 25), middleage and old (over 50).}
\item{training}{Whether an individual received a professional training.}
\item{university}{Whether an individual has a university degree or comparable.}
}
}

\source{Socio Economic Panel (SOEP), s. \url{www.diw.de/deutsch/sop}.}
\references{Kauermann G. and Khomski P. (2006). Additiv two way hazards model with varying coefficients, in press.
}
\keyword{datasets}

