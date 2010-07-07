\name{addManual}
\alias{addManual}
\alias{delManual}
\title{ Add or delete stations manually to/from a given design }
\description{ Add or delete stations manually to/from a given design }
\usage{
addManual(candidates, observations, nDiff)
delManual(candidates, observations, nDiff)
}
\arguments{
\item{observations}{ Data frame with x,y coordinates of observations and possible covariates }
\item{candidates}{ When action is "add", candidates is data frame with x,y coordinates }
\item{nDiff}{number of points to add or delete}
}

\value{
Data frame with added or deleted observations optimized locations
}

\references{
\url{http://www.intamap.org/}
}
\author{ Olivier Baume }
\keyword{}
