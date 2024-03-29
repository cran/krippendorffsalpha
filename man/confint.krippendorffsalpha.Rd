% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/krippendorff.R
\name{confint.krippendorffsalpha}
\alias{confint.krippendorffsalpha}
\title{Compute a confidence interval for Krippendorff's Alpha.}
\usage{
\method{confint}{krippendorffsalpha}(object, parm = "alpha", level = 0.95, ...)
}
\arguments{
\item{object}{an object of class \code{"krippendorffsalpha"}, the result of a call to \code{\link{krippendorffs.alpha}}.}

\item{parm}{always ignored since there is only one parameter.}

\item{level}{the desired confidence level for the interval. The default is 0.95.}

\item{\dots}{additional arguments. These are passed to \code{\link{quantile}}.}
}
\value{
A vector with entries giving lower and upper confidence limits. These will be labelled as (1 - level) / 2 and 1 - (1 - level) / 2.
}
\description{
Compute a confidence interval for Krippendorff's Alpha.
}
\details{
This function computes a confidence interval for alpha, assuming that \code{\link{krippendorffs.alpha}} was called with \code{confint = TRUE}.

For \code{method = "analytical"}, a jackknife-based interval is computed. For smaller samples the jackknife interval offers a very substantial improvement over the bootstrap interval, the latter of which offers quite poor coverage. For larger samples \code{method = "customary"} can safely be used, in which case a bootstrap interval is provided. For sufficiently large datasets the two intervals will be nearly equal, but the bootstrap approach is preferred owing to its much faster execution speed.
}
\examples{
# Fit a subset of the cartilage data, using the customary methodology.
# Compute bootstrap confidence intervals using a bootstrap sample size
# of 1,000. Report the estimate of alpha, and produce a 99\% interval.

data(cartilage)
cartilage = as.matrix(cartilage[1:100, ])
fit.cart = krippendorffs.alpha(cartilage, level = "ratio", method = "customary", confint = TRUE,
                               control = list(bootit = 1000, parallel = FALSE))
fit.cart$alpha.hat
confint(fit.cart, level = 0.99)
}
\references{
Nissi, M. J., Mortazavi, S., Hughes, J., Morgan, P., and Ellermann, J. (2015). T2* relaxation time of acetabular and femoral cartilage with and without intra-articular Gd-DTPA2 in patients with femoroacetabular impingement. \emph{American Journal of Roentgenology}, \bold{204}(6), W695.
}
\seealso{
\code{\link{krippendorffs.alpha}}
}
