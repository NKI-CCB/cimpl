\name{doCimplAnalysis}
\alias{doCimplAnalysis}
\title{
Do "cimpl" analysis
}
\description{
Calculate a \linkS4class{cimplAnalysis} object using insertion data.
}
\usage{
doCimplAnalysis(data, scales = (1:15) * 10000, n_iterations = 100, sample_period = 0.1, max_sample_points = 2^19, chromosomes = unique(data$chr), BSgenome, system = c("MMTV", "MuLV", "SB", "PB"), specificity.pattern, lhc.method = c("none", "exclude"), verbose = TRUE)
}

\arguments{
  \item{data}{A data.frame which contains the insertion locations by row. Required is a 'location' column and a 'chr' (chromosome) column, allowing only 'chr1', ..., 'chrY' values. If \code{lhc.method = 'exclude'} a 'contig_depth' column is required. Optional columns are 'tumorID' and 'sampleID' which are used in \code{plot(..., type='kse+')}. Not required columns are ignored, but are maintained in \code{\linkS4class{cimplObject}@data}.}
  \item{scales}{The range of scales used for the analysis.}
  \item{n_iterations}{The number of random permutations used for calculation of the null-distribution of peak heights.}
  \item{sample_period}{The sample period as a percentage of the scale. Used to calculate the sample points at which the densities are estimated.}
  \item{max_sample_points}{The maximum number of sample points.}
  \item{chromosomes}{The chromosomes used for this analysis. Formatted as "chr1", ..., "chrY".}
  \item{BSgenome}{The genome used in the experiment. For example: 'Mmusculus' from the BSgenome.Mmusculus.UCSC.mm9 package. NB: The 'data' refers to this genome!}
  \item{system}{The system used in the experiment. Possible values are the virus-based systems 'MMTV' (Mouse Mammary Tumor Virus) and 'MuLV' (Murine Leukemia Virus) and the transposon based systems 'SB' (Sleeping Beauty) and 'PB' (piggyBac).}
  \item{specificity.pattern}{The specificity pattern at which the insertion locations are restricted. If no system is set, the specificity pattern can be user-defined.}
  \item{lhc.method}{The local hopping correction method. Possible values are 'none' or 'exclude'. When the distance between two neighboring insertions is less than three kernel widths, the insertion with the smallest 'contig_depth' is the hopped insertion. If \code{lhc.method == 'exclude'}, hopped insertions are excluded from the analysis.}
  \item{verbose}{Echo progress indicators?}
  \item{cores}{Number of cores to use when calculating the null model.}
}
\details{

If \code{specificity.pattern != ''}, background densities are calculated using Gaussian kernel density estimation for each scale. (Note: else the background density is uniform)

Then, for each chromosome and each scale:

1. The data is randomized, i.e each insertion (if \code{lhc.method=='exclude'} hops are excluded) is given a random location on the chromosome. If \code{specificity.pattern != ''} possible locations are restricted accordingly.

2. Gaussian kernel density estimation is applied on the randomized data and the peaks are located. For each peak, the density (or peak height) and the background density at the peak location are added to the (joint) null-distribution.

3. Steps 1 and 2 are repeated \code{n_iteration} times.

4. Gaussian kernel density estimation is applied on the observed (\code{data}). Peaks are located, and given a significance value using the generated null-distribution.

}
\value{
	The results are placed in a \linkS4class{cimplAnalysis} object.
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
}
\seealso{
\code{\linkS4class{cimplAnalysis}}
}

