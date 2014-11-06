\name{getInsertions}
\alias{getInsertions}
\title{
	Get Insertions
}
\description{
	Returns the insertions from a given cimplAnalysis object.
}
\usage{
	getInsertions(cimplAnalysis, chr, scale, bpLim)
}
\arguments{
	\item{cimplAnalysis}{The \linkS4class{cimplAnalysis} object.}
	\item{chr}{The chromosome.}
	\item{scale}{The scale.}
	\item{bpLim}{The region on the chromosome: starting at \code{bpLim[1]} and ending at \code{bpLim[2]}.}
}
\value{
	A data.frame of insertions that meet the requirements set by the parameters. The \code{data} object provided in \link{doCimplAnalysis} is subsetted by row, i.e. the resulting data.frame contains the same columns. Additionally, if a local hopping correction method was applied, a logical column \code{isHop} is appended.
}
\author{
	Jelle ten Hoeve
}
\examples{
\dontrun{
library(BSgenome.Mmusculus.UCSC.mm9)
data(colorectal)

# do 'cimpl' analysis
sampleCa <- doCimplAnalysis(colorectal, scales = c(100e3), chromosomes = c('chr19'), n_iterations = 100, BSgenome = Mmusculus, system = 'SB', lhc.method='none')
}
data(sampleCa) # load sample data to reduce package build time

insertions <- getInsertions(sampleCa, chr='chr18', scale=50e3, bpLim=c(1, 10e6))

# show 10 insertions
insertions[1:10,]
}