\name{simdata}
\docType{data}
\alias{simdata}
\title{Simulated data}
\description{
Data simulation using the LEMMA model.  The true parameter values are:
\eqn{\alpha=5, \beta=1/12, \tau=0, \psi=-1, \sigma^2_\psi=0.5, p_1=p_2=0.25}.
Both groups have 6 replicates (`subjects'), and 2,000 `genes'.
}
\usage{simdata}
\format{A data frame containing 2,000 rows and 14 columns.  The first column contains
the gene IDs, the second column contains the gene names.  Columns 3-8 contain the
reponses for the 6 control samples, and columns 9-14 contain the
reponses for the 6 treatment samples.}
\keyword{datasets}
