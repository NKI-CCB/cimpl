\name{cimplAnalysis-class}
\docType{class}
\alias{cimplAnalysis-class}
\alias{show,cimplAnalysis-method}

\title{Class "cimplAnalysis"}
\description{Contains the results of a "cimpl" analsysis.}
\section{Objects from the Class}{
Objects can be created by calls of the method \link{doCimplAnalysis}.
}
\section{Slots}{
	 \describe{
    \item{\code{cimplObjects}:}{A two dimensional list (a list of lists) containing \linkS4class{cimplObject} objects. To obtain a \linkS4class{cimplObject} object for a specific chromosome and scale, use \code{...@cimplObjects[[chromosome.idx]][[scale.idx]]}.}
    \item{\code{chromosomes}:}{The chromosomes used for this analysis. Formatted as "chr1", ..., "chrY".}
    \item{\code{scales}:}{The scale used for this analysis.}
    \item{\code{system}:}{The system used in the experiment. Possible values are the virus-based systems 'MMTV' (Mouse Mammary Tumor Virus) and 'MuLV' (Murine Leukemia Virus) and the transposon based systems 'SB' (Sleeping Beauty) and 'PB' (piggyBac).}
    \item{\code{specificity.pattern}:}{The specificity pattern at which the insertion locations are restricted. If no system is set, the specificity pattern can be user-defined.}
    \item{\code{lhc.method}:}{The local hopping correction method. Possible values are 'none' or 'exclude'.}
    \item{\code{n_iterations}:}{The number of random permutations used for calculation of the null-distribution of peak heights.}
    \item{\code{sample_period}:}{The sample period as a percentage of the scale. Used to calculate the sample points at which the densities are estimated.}
    \item{\code{max_sample_points}:}{The maximum number of sample points.}
    \item{\code{data.name}:}{The data used for this analysis.}
    \item{\code{call}:}{The original 'call' of which this \linkS4class{cimplAnalysis} object is the result.}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "cimplAnalysis", y = "missing")}: \code{\link{plot}}}
    \item{show}{\code{signature(object = "cimplAnalysis")}: \code{\link{show}}}
	 }
}

\seealso{
\code{\link{doCimplAnalysis}}, \code{\linkS4class{cimplObject}}
}
\examples{
showClass("cimplAnalysis")
}
\keyword{classes}
