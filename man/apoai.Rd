\name{apoai}
\docType{data}
\alias{apoai}
\title{The APO-AI data set (Callow et al., 2000) }
\description{
The ApoA1 experiment (Callow et al., 2000) used gene targeting in embryonic
stem cells to produce mice lacking apolipprotein A-1, a gene known to play a
critical role in high density lipoprotein (HDL) cholesterol levels. Originally,
5,600 expressed sequence tags (EST) were selected. For this data set, we used
the data and normalization method provided with the LIMMA R package
(Smyth, 2005), which consists of 5,548 ESTs, from 8 control (wild type "black
six") mice and 8 "knockout" (lacking ApoA1) mice. Common reference RNA
was obtained by pooling RNA from the control mice, and was used to perform
expression profiling for all 16 mice.
The response of interest  is the log2 fluorescence ratio (with respect
to the common reference).
}
\usage{apoai}
\format{A data frame containing 5,548 rows and 18 columns.  The first column contains
the gene IDs, the second column contains the gene names.  Columns 3-11 contain the
reponses for the 8 mice in the control group, and columns 12-18 contain the
reponses for the 8 mice in the treatment group.}
\references{
Callow, M. J., Dudoit, S., Gong, E. L., Speed, T. P. and Rubin,
E. M. (2000). Microarray expression profiling identifies genes with altered
expression in HDL-deficient mice. Genome Res, 10 2022-9.
}
\keyword{datasets}
