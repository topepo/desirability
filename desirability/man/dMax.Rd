\name{dMax}
\alias{dMax}
\alias{dMin}
\alias{dTarget}
\alias{dArb}
\alias{dBox}
\alias{dCategorical}
\alias{dOverall}

\alias{dMax.default}
\alias{dMin.default}
\alias{dTarget.default}
\alias{dArb.default}
\alias{dBox.default}
\alias{dCategorical.default}
\alias{dOverall.default}

\title{Desirability Functions}
\description{
Functions implementing multivariate optimization and ranking using the desirability function approach described in Derringer and Suich (1980)
}
\usage{
\method{dMax}{default}(low, high, scale = 1, tol = NULL, ...)
\method{dMin}{default}(low, high, scale = 1, tol = NULL, ...)
\method{dTarget}{default}(low, target, high, lowScale = 1, highScale = 1, tol = NULL, ...) 
\method{dArb}{default}(x, d, tol = NULL, ...)
\method{dBox}{default}(low, high, tol = NULL, ...)
\method{dCategorical}{default}(values, tol = NULL, ...)
\method{dOverall}{default}(...)

}
\arguments{
  \item{low}{a constant to define the desirability function for \code{dMax}, \code{dMin}, \code{dTarget} and \code{dBox}}
  \item{high}{a constant to define the desirability function for \code{dMax}, \code{dMin}, \code{dTarget} and \code{dBox}}
  \item{target}{a constant to define the desirability function for \code{dMax}, \code{dMin}, \code{dTarget} and \code{dBox}}
  \item{scale}{the scaling factor for \code{dMax} and \code{dMin}. Values less than one make the criteria more difficult to satisfy while values greater than one make it easier.} 
  \item{lowScale}{the scaling factor for \code{dTarget}. This bends the curve between the points \code{low} and \code{target}. Values less than one make the criteria more difficult to satisfy while values greater than one make it easier. }
  \item{highScale}{ the scaling factor for \code{dTarget}. This bends the curve between the points \code{high} and \code{target}. Values less than one make the criteria more difficult to satisfy while values greater than one make it easier.}  
  \item{x}{a set of input values}
  \item{d}{a set of desirabilites between zero and one (inclusive) that match the length of \code{x}}    
  \item{values}{a named numeric vector of possible values}
  \item{tol}{an optional tolerance for zero desirability. When this is non-null, zero desirabilites are replaced with this value}
  \item{\ldots}{For \code{dOverall}, this is one or more desirability objects. For the other methods, this argument is not currently used}
}
\details{
The functions  \code{dMax}, \code{dMin}, \code{dTarget} ande \code{dOverall} are the basic equations used by Derringer and Suich (1980). \code{dBox} is a simple step funciton between two points. \code{dArb} can be used to create other shapes that do not fall into the other funcional forms. See the package vignette or the references for more details
}
\value{
a list. Common values are:

  \item{tol }{the value specified by the \code{tol} argument}
  \item{call }{the origianl function call}

}
\references{Derringer, G. and Suich, R. (1980),
Simultaneous Optimization of Several Response Variables.
\emph{Journal of Quality Technology} \bold{ 12}, 214--219.}
\author{Max Kuhn}

\seealso{\code{\link{predict.dMax}}}
\examples{
dMax.default(1,3)
dMax(1,3)
}
\keyword{utilities}

