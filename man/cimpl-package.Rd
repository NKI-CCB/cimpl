\name{cimpl-package}
\alias{cimpl-package}
\alias{cimpl}
\docType{package}
\title{
Common Insertion site Mapping PLatform
}
\description{
An analysis package for multi sample insertional mutagenesis data (including viral- and transposon-based systems) using Gaussian kernel convolution to identify common insertion sites.
}
\author{
Jelle ten Hoeve and Jeroen de Ridder
}
\references{
De Ridder J, Uren A, Kool J, Reinders M, and Wessels L. Detecting statistically significant common insertion sites in retroviral insertional mutagenesis screens. PLoS Comput Biol, 2:e166, 2006.
}

\keyword{ package }
\seealso{ 
\code{\link{doCimplAnalysis}}, \code{\link{plot}}
}
\examples{
\dontrun{
library(BSgenome.Mmusculus.UCSC.mm9)
data(colorectal)

# do 'cimpl' analysis
sampleCa <- doCimplAnalysis(colorectal, scales = c(100e3), chromosomes = c('chr19'), n_iterations = 100, BSgenome = Mmusculus, system = 'SB', lhc.method='none')
}
data(sampleCa) # load sample data to reduce package build time

# make some plots
plot(sampleCa, type='kse', interactive=FALSE)
plot(sampleCa, type='scale.space', interactive=FALSE)
plot(sampleCa, type='null.cdf', interactive=FALSE)

genes <- getEnsemblGenes(sampleCa)

# retrieve CISs
ciss <- getCISs(sampleCa, genes=genes, alpha=0.05, mul.test=TRUE)  # NB: set significance level (alpha) and multiple testing correction
write.csv(ciss, file='ciss.csv')

# export result to html
export.html(sampleCa, genes=genes, alpha=0.05, mul.test=TRUE)

# make a matrix linking insertions to CISs
cisMatrix <- getCISMatrix(sampleCa, ciss)
write.csv(cisMatrix, file='cisMatrix.csv', row.names=FALSE)

# export kse and bg_density to .wig file
export.wig(sampleCa, file='kse.wig')

# export CISs to .bed file
export.bed(ciss, file='ciss.bed')

}
