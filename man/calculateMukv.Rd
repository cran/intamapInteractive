\name{calculateMukv}
\alias{calculateMukv}
\title{MUKV computation}
\description{
Computes MUKV for given geostatistical parameters
}
\usage{
calculateMukv(observations, predGrid, model, formulaString, ...)
}
\arguments{
\item{observations}{ \code{\link[sp]{SpatialPoints}} or  \code{\link[sp]{SpatialPointsDataFrame}} 
     with observation locations and possible covariates }
\item{predGrid}{ \code{\link[sp:Spatial-class]{Spatial}} object with coordinates of 
     prediction locations (usually \code{\link[sp:SpatialPixels]{SpatialGrid}} 
       or \code{\link[sp:SpatialPixelsDataFrame]{SpatialGridDataFrame}})}
\item{model}{ Variogram model:object of class \code{variogramModel}, of the form
                created by \code{\link[gstat:vgm]{vgm}}}
\item{formulaString}{formula that defines the dependent variable as a linear model 
      of independent variables; suppose the dependent variable has name \code{z}, 
      for ordinary and simple kriging use the formula \code{z~1}; 
      for universal kriging, suppose \code{z} is linearly dependent on 
      \code{x} and \code{y}, use the formula \code{z~x+y}. The formulaString defaults
      to \code{"value~1"} if \code{value} is a part of the data set. 
      If not, the first column of the data set is used.}
\item{...}{ other arguments to be passed on at lower level functions }
}

\details{
This function computes kriging on the \code{predGrid} with 
\code{\link[gstat:krige]{krige}} function, and averages the kriging variance (MUKV). With covariates, 
the function takes a universal kriging model into account.
}

\value{
MUKV value
}

\references{
\url{http://www.intamap.org/}
}
\author{ S.J. Melles and O. Baume }
\keyword{spatial}
