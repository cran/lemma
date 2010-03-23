\name{ecoli}
\docType{data}
\alias{ecoli}
\title{The Escherichia coli data (Hung et al., 2002) }
\description{
The Leucine-responsive regulatory protein (Lrp) data analyzed by Hung et al.
consists of 2,758 genes, with 4 replicates in each treatment group.
The experimental design for the Affymetrix GeneChip experiments is as follows:
four GeneChips were hybridized with biotin-labeled RNA pools 1-3, 4-6, 7-9,
and 10-12 prepared from lrp+ cells, and four GeneChips were hybridized with
biotin-labeled RNA pools 1-3, 4-6, 7-9, and 10-12 prepared from lrp- cells,
respectively. One average difference measurement for each gene probe set on
each GeneChip was obtained for subsequent data processing and analysis. 
The normalized data was downloaded from
\cr http://www.jbc.org/cgi/content/full/277/43/40309/DC1 and was further
transformed by taking the log2 of the responses.
}
\usage{ecoli}
\format{A data frame containing 2,758 rows and 10 columns.  The first column
contains the gene IDs, the second column contains the gene names.
Columns 3-6 contain the reponses for the 4 control samples, and columns 7-10
contain the reponses for the 4 treatment samples.}
\references{
Hung, S. P., Baldi, P. & Hatfield, G. W. (2002). Global gene expression
profiling in Escherichia coli K12. The effects of leucine-responsive
regulatory protein.
J Biol Chem 277, 40309-40323
}
\keyword{datasets}
