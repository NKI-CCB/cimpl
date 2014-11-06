\name{colorectal}
\alias{colorectal}
\docType{data}
\title{
	Colorectal Cancer Sleeping Beauty dataset
}
\description{
	Data retrieved from 135 (cohort 1) and 93 (cohort 2) tumors induced using a Sleeping Beauty transposon-based system.
}
\usage{data(colorectal)}
\format{
  A data frame with 16690 non-redundant mapped insertions.
  \describe{
    \item{\code{sampleID}}{The mouse identifier.}
    \item{\code{tissueType}}{Section of GI tract: Col = colon; Duo = duodenum, Jej = Jejunum, Ile = ileum, Cec = cecum, Int = small intestine (exact location not specified)}
    \item{\code{tissueID}}{When more than one tumor was taken from the same section of intestine, tumors were labeled consecutively, (i.e. A, B, C, etc.).}
    \item{\code{chr}}{Chromosome}
    \item{\code{location}}{Chromosomal address corresponds to NCBI Mouse Build 37.}
    \item{\code{associatedGene}}{If insertion contributed to a CIS, candidate mouse gene is listed in this column.}
    \item{\code{cohort}}{Dataset 1 = mice from first cohort, dataset 2 = mice from second cohort.}
    \item{\code{contig_depth}}{The read count. NB: all set 1 since read count is not published.}
  }
}
\source{
A transposon-based genetic screen in mice identifies genes altered in colorectal cancer.
Starr TK, Allaei R, Silverstein KA, Staggs RA, Sarver AL, Bergemann TL, Gupta M, O'Sullivan MG, Matise I, Dupuy AJ, Collier LS, Powers S, Oberg AL, Asmann YW, Thibodeau SN, Tessarollo L, Copeland NG, Jenkins NA, Cormier RT, Largaespada DA.
Science. 2009 Mar 27;323(5922):1747-50. Epub 2009 Feb 26.
}
\examples{
data(colorectal)
str(colorectal)
}
\keyword{datasets}
