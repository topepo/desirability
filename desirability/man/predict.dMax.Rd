\name{predict.dMax}
\alias{predict.dMax}
\alias{predict.dMin}
\alias{predict.dTarget}
\alias{predict.dArb}
\alias{predict.dBox}
\alias{predict.dOverall}

\title{Predict method for desirability functions}
\description{
Predicted values based on desirability objects
}
\usage{
\method{predict}{dMax}(object, newdata = NA, ...)
\method{predict}{dMin}(object, newdata = NA, ...)
\method{predict}{dTarget}(object, newdata = NA, ...)
\method{predict}{dArb}(object, newdata = NA, ...)
\method{predict}{dBox}(object, newdata = NA, ...)
\method{predict}{dOverall}(object, newdata = matrix(NA, ncol = length(object$d)), all = FALSE, ...)

}
\arguments{
  \item{object}{a object of class: \code{dMax}, \code{dMin},
    \code{dTarget}, \code{dArb}, \code{dBox} or \code{dOverall}}
  \item{newdata}{values of the response for predicting desirability}
  \item{all}{a logical (for \code{predict.dOverall} only); should the
    individual desirabilities also be returned?}
  \item{\dots}{no currently used}
}
\details{
The responses are translated into desirability units.
}
\value{
a vector, unless \code{predict.dOverall} is used with \code{all=TRUE},
in which case a matrix is returned.
}
\references{Derringer, G. and Suich, R. (1980),
Simultaneous Optimization of Several Response Variables.
\emph{Journal of Quality Technology} \bold{ 12}, 214--219.}
\author{Max Kuhn}

\seealso{\code{\link{dMax}}}
\examples{
d1 <- dMin(1,3)
d2 <- dTarget(1, 2, 3)
dBoth <- dOverall(d1, d2)

outcomes <- cbind(
   seq(0, 4, length = 10),
   seq(0.5, 4.5, length = 10))
   
  
predict(d1, outcomes[,2])
predict(d2, outcomes[,2])
predict(dBoth, outcomes)
predict(dBoth, outcomes, all = TRUE)
}
\keyword{utilities}

