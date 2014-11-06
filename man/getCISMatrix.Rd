\name{getCISMatrix}
\alias{getCISMatrix}
\title{
	Get CIS matrix
}
\description{
	Returns a CIS matrix, which indicates the CISs associated with insertions.
}
\usage{
getCISMatrix(cimplAnalysis, ciss)
}
\arguments{
  \item{cimplAnalysis}{A \linkS4class{cimplAnalysis} object.}
  \item{ciss}{A \code{data.frame} containing the CISs (\link{getCISs}).}
}
\value{
	A data.frame where rows are the insertions. The first columns contain information about the insertions, these are exactly the same columns as provided in \code{data} (\link{doCimplAnalysis}. After that, for each scale a column is added containing CIS ids, which can be trace back in the provided \code{ciss} data.frame.).  
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

genes <- getEnsemblGenes(sampleCa)

# retrieve CISs
ciss <- getCISs(sampleCa, genes=genes, alpha=0.05, mul.test=TRUE)  # NB: set significance level (alpha) and multiple testing correction

# make a matrix linking insertions to CISs
cisMatrix <- getCISMatrix(sampleCa, ciss)
write.csv(cisMatrix, file='cisMatrix.csv', row.names=FALSE)
}