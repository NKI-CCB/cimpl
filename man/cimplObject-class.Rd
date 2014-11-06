\name{cimplObject-class}
\docType{class}
\alias{cimplObject-class}
\alias{show,cimplObject-method}

\title{Class "cimplObject"}
\description{Represents the partial "cimpl" analysis of a single chromosome and a single scale. }
\section{Objects from the Class}{
Objects are created in \link{doCimplAnalysis} and can be retrieved form a \linkS4class{cimplAnalysis} object using \code{...@cimplObjects[[chromosome.idx]][[scale.idx]]}.
}
\section{Slots}{
	 \describe{
    \item{\code{data}:}{A data.frame of the insertions on \code{chromosome}. All the columns provided in \link{doCimplAnalysis} are included and if \code{lhc.method != 'none'} a logical column 'isHop' is provided.}
    \item{\code{kse}:}{A list containing a \code{x} and \code{y} vector, representing the Gaussian kernel density estimation of the observed \code{data} using \code{scale}.}
    \item{\code{peaks}:}{A list containing a \code{x}, \code{y} and \code{p-value} vector representing the peak coordinates and significances.}
    \item{\code{null_peaks}:}{A list containing a \code{x} and \code{y} vector, representing all the peaks in the null-distribution.}
    \item{\code{null_cdf}:}{A list containing a \code{x}, \code{y} and \code{z} vector, representing the joint peak-height/background-density distribution, conditioned on background-density.}
    \item{\code{bg_density}:}{A list containing a \code{x} and \code{y} vector, representing the background density.}
    \item{\code{n_peaks}:}{The number of peaks.}
    \item{\code{scale}:}{The scale of this partial analysis.}
    \item{\code{chromosome}:}{The chromosome of this partial analysis.}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "cimplObject", y = "missing")}: \code{\link{plot}} }
    \item{show}{\code{signature(object = "cimplObject")}: \code{\link{show}}}
	 }
}
\seealso{
\code{\linkS4class{cimplAnalysis}}, \code{\link{doCimplAnalysis}}
}
\examples{
showClass("cimplObject")
}
\keyword{classes}
