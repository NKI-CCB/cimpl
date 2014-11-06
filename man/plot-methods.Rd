\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,cimplAnalysis,missing-method}
\alias{plot,cimplObject,missing-method}

\title{Plot "cimplAnalysis" or "cimplObject"}
\description{
Methods for function \code{plot} in package `cimpl'
}
\usage{
\S4method{plot}{cimplAnalysis,missing}(x, y, type=c('kse', 'null.cdf', 'scale.space'), chr=x@chromosomes[1], scale=x@scales[1], alpha=0.05, mul.test=TRUE, bpLim, plot.tumor.densities=FALSE, interactive=TRUE, ...)
\S4method{plot}{cimplObject,missing}(x, y, type=c('kse', 'null.cdf'), alpha = 0.05, n_tests=1, bpLim, plot.tumor.densities=FALSE, interactive=TRUE, ...)
}
\arguments{
  \item{x}{An object of class cimplAnalysis or cimplObject.}
  \item{y}{An object of class missing (do not provide!).}
  \item{type}{The plot type: 'kse' for the kernel smoothed estimate of the data,.}
  \item{chr}{The chromosome ('chr1', ..., 'chrY').}
  \item{scale}{The scale.}
  \item{alpha}{The significance level, plotted as solid line.}
  \item{mul.test}{Correct for multiple testing? A dotted red line is added if \code{TRUE}.}
  \item{n_tests}{The total number of tests. Automatically calculated when plotting a cimplAnalysis and \code{mul.test == TRUE}.}
  \item{bpLim}{The range of the x-axis.}
  \item{plot.tumor.densities}{If \code{TRUE}, the contribution of each tumor to the kse (if 'sampleID' or 'tumorID' is supplied in \code{data} (\link{doCimplAnalysis})) is visualized in the plot.}
  \item{interactive}{If \code{TRUE}, the user can navigate and zoom using the mouse. Set to \code{FALSE} if you plot to a file.}
  \item{...}{All other arguments are passed to \code{plot}.}
}
\keyword{methods}
\examples{
\dontrun{
library(BSgenome.Mmusculus.UCSC.mm9)
data(colorectal)

# do 'cimpl' analysis
sampleCa <- doCimplAnalysis(colorectal, scales = c(100e3), chromosomes = c('chr19'), n_iterations = 100, BSgenome = Mmusculus, system = 'SB', lhc.method='none')
}
data(sampleCa) # load sample data to reduce package build time

plot(sampleCa, type='kse', interactive=FALSE)
plot(sampleCa, type='scale.space', interactive=FALSE)
plot(sampleCa, type='null.cdf', interactive=FALSE)

}