\name{predict.dMax}
\alias{predict.dMax}
\alias{predict.dMin}
\alias{predict.dTarget}
\alias{predict.dArb}
\alias{predict.dBox}
\alias{predict.dCategorical}
\alias{predict.dOverall}

\title{Predict method for desirability functions}
\description{
Predicted values based on desirability objects
}
\usage{
\method{predict}{dMax}(object, newdata = NA, missing = object$missing, ...)
\method{predict}{dMin}(object, newdata = NA, missing = object$missing, ...)
\method{predict}{dTarget}(object, newdata = NA, missing = object$missing, ...)
\method{predict}{dArb}(object, newdata = NA, missing = object$missing, ...)
\method{predict}{dBox}(object, newdata = NA, missing = object$missing, ...)
\method{predict}{dCategorical}(object, newdata = NA, missing = object$missing, ...)
\method{predict}{dOverall}(object, newdata = data.frame(NA, ncol = length(object$d)), all = FALSE, ...)

}
\arguments{
  \item{object}{a object of class: \code{dMax}, \code{dMin},
    \code{dTarget}, \code{dArb}, \code{dBox} or \code{dOverall}}
  \item{newdata}{values of the response for predicting desirability}
  \item{all}{a logical (for \code{predict.dOverall} only); should the
    individual desirabilities also be returned?}
  \item{missing}{a number between 0 and 1 for missing values (the internally estimated value is used by default)}
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
d3 <- dCategorical(c("a" = .1, "b" = .7))
dAll <- dOverall(d1, d2, d3)

outcomes <- data.frame(seq(0, 4, length = 10),
                       seq(0.5, 4.5, length = 10),
                       sample(letters[1:2], 10, replace = TRUE))
names(outcomes) <- c("x1", "x1", "x3")   
  
predict(d1, outcomes[,2])
predict(d2, outcomes[,2])
predict(d3, outcomes[,3])
predict(dAll, outcomes)
predict(dAll, outcomes, all = TRUE)
}
\keyword{utilities}

