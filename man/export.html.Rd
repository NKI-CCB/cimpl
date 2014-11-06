\name{export.html}
\alias{export.html}
\alias{export.bed}
\alias{export.wig}
\title{
	Export the results
}
\description{
	Generate HTML files to inspect the results within your browser.
	Generate WIG files to include the densities in a genome browser like Ensembl Genome Browser.
	Generate BED files to include the CISs in a genome browser.
}
\usage{
export.html(cimplAnalysis, chromosomes = cimplAnalysis@chromosomes, scales = cimplAnalysis@scales, genes, alpha = 0.05, mul.test = TRUE, plot.tumor.densities = FALSE, order.by = c("scale", "n_insertions", "p_value"), decreasing = c(FALSE, TRUE, FALSE), dir = ".", verbose = TRUE)
export.wig(cimplAnalysis, chromosomes = cimplAnalysis@chromosomes, scales = cimplAnalysis@scales, file)
export.bed(ciss, file)
}
\arguments{
  \item{cimplAnalysis}{A \linkS4class{cimplAnalysis} object.}
  \item{chromosomes}{By default all available chromosomes.}
  \item{scales}{By default all available scales.}
  \item{genes}{If missing, genes are retrieved for BioMart.}
  \item{alpha}{The significance threshold used for obtaining the CISs.}
  \item{mul.test}{Should multiple testing be applied?}
  \item{plot.tumor.densities}{Should the tumor densities by included in the kse plots?}
  \item{order.by}{The ordering of the CISs.}
  \item{decreasing}{The direction of the ordering.}
  \item{file}{The output file.}
  \item{dir}{The output directory.}
  \item{verbose}{Echo progress indicators?}
  \item{ciss}{A \code{data.frame} containing the CISs (\link{getCISs}).}
}
\value{
	HTML files and image files are placed in the output directory (\code{dir}). Open index.html with your browser.
	WIG files and BED files are generated using the \code{file} that is supplied.
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

# export result to html
export.html(sampleCa, genes=genes, alpha=0.05, mul.test=TRUE)

# export kse and bg_density to .wig file
export.wig(sampleCa, file='kse.wig')

# export CISs to .bed file
ciss <- getCISs(sampleCa, genes=genes, alpha=0.05, mul.test=TRUE)  # NB: set significance level (alpha) and multiple testing correction
export.bed(ciss, file='ciss.bed')
}

