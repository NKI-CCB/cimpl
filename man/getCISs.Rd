\name{getCISs}
\alias{getCISs}
\title{
	Get Common Insertion Sites
}
\description{
	A CIS is a genomic region for which the 'kse' is above the significance threshold. For each CIS found, an algorithm is ran associate genes to it.
}
\usage{
	getCISs(cimplAnalysis, alpha = 0.05, chromosomes = cimplAnalysis@chromosomes, scales = cimplAnalysis@scales, mul.test = TRUE, mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl"), genes = getEnsemblGenes(cimplAnalysis, mart), order.by = c("p_value", "n_insertions"), decreasing = c(FALSE, TRUE))
}
\arguments{
	\item{cimplAnalysis}{A \linkS4class{cimplAnalysis} object.}
	\item{alpha}{The significance level.}
	\item{chromosomes}{Only CISs on these \code{chromosomes} are returned.}
	\item{scales}{Only CISs for these \code{scales} are returned.}
	\item{mul.test}{A logical indicating if a mulitple testing (Bonferroni) correction has to be applied.}
	\item{mart}{The mart (see biomaRt package) from which the genes are retrieved if \code{genes} are missing.}
	\item{genes}{A data.frame containing genes. For details: \link{getEnsemblGenes}.}
	\item{order.by}{A vector of column names (see Value section) by which the resulting CIS list is ordered.}
	\item{decreasing}{A logical vector of same length as \code{order.by}.}
}
\value{
	A data.frame of CISs containing the following columns. If \code{collapse == FALSE}, CISs associated with more than one gene will span multiple rows and the 'other_genes' column is omitted.
	\item{external_gene_id}{The gene closest to the (highest if more than 1) peak of the CIS and overlapping with the CIS. If no gene is found, this field is blank ('').}
	\item{chromosome}{The chromosome of the CIS.}
	\item{peak_location}{The location of the (highest if more than 1) peak within the CIS.}
	\item{start}{The start address of the CIS.}
	\item{end_pos}{The end address of the CIS.}
	\item{width}{The CIS width in base pairs.}
	\item{n_insertions}{The number of insertions found within the CIS.}
	\item{other_genes}{Other genes found to overlap with the CIS. If \code{collapse == FALSE}, this column is omitted.}
	\item{p_value}{The significance value of the peak with largest peak height within the significant region. The data.frame is ordered by increasing p-value.}
	\item{scale}{The scale at which the CIS is found.}
}
\author{
	Jelle ten Hoeve
}
\seealso{
	\link{getInsertions}
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

write.csv(ciss, file='ciss.csv')
}
