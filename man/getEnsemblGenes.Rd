\name{getEnsemblGenes}
\alias{getEnsemblGenes}
\title{
	Get Ensembl Genes
}
\description{
	Retrieves Ensembl genes using the biomaRt package.
}
\usage{
getEnsemblGenes(cimplAnalysis, mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl"))
}
\arguments{
	\item{cimplAnalysis}{ The \linkS4class{cimplAnalysis} object.}
	\item{mart}{A mart object (see biomaRt package).}
}
\value{
	A data.frame of genes, containing the 'external_gene_id' , 'chromosome_name', 'start_position' and 'end_position' attributes.
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

# show 10 genes
genes[1:10,]
}
