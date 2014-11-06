\name{sampleCa}
\alias{sampleCa}
\docType{data}
\title{
	Sample CimplAnalysis object
}
\description{
	Used in the examples of this package.
}
\usage{data(sampleCa)}
\examples{
\dontrun{
library(BSgenome.Mmusculus.UCSC.mm9)
data(colorectal)

# do 'cimpl' analysis
sampleCa <- doCimplAnalysis(colorectal, scales = c(100e3), chromosomes = c('chr19'), n_iterations = 100, BSgenome = Mmusculus, system = 'SB', lhc.method='none')
}

data(sampleCa)
str(sampleCa)
}
\keyword{datasets}
