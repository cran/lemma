\name{colon}
\docType{data}
\alias{colon}
\title{The colon cancer data set (Alon et al., 1999) }
\description{
The data analyzed by Alon et al. (1999) consists of 2000 ESTs in 40 tumor
and 22 normal colon tissue samples. Of the 40 patients involved in the
study, 22 supplied both tumor and normal tissue samples.
The normalized data was downloaded from \cr http://microarray.princeton.edu/oncology/affydata/index.html and was further transformed by taking the log2 of the expression data.
}
\usage{colon}
\format{A data frame containing 2,000 rows and 64 columns.  The first column contains
the gene IDs, the second column contains the gene names.  Columns 3-24 contain the
reponses for the 22 control samples, and columns 25-64 contain the
reponses for the 40 treatment samples.}
\references{
Alon, Barkai, Notterman, Gish, Ybarra, Mack and Levine (1999).
Broad patterns of gene expression revealed by clustering analysis of tumor
and normal colon tissues probed by oligonucleotide arrays. Proceedings of
the National Academy of Sciences, 96 6745-6750.
}
\keyword{datasets}
