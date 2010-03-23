\name{lemmaPlots}
\alias{lemmaPlots}
\docType{methods}
\title{LEMMA - diagnostics plots}
\description{
This program produces five diagnostics plots in two files (in postscript and
pdf formats).  The file LEMMAplots1 contains four plots:
(a) the fitted distribution for the histogram of \eqn{d_g}, (b) the
fitted distribution for the histogram of \eqn{m_g}, (c) and (d) are
`volcano plots' that show the distribution of the Benjamini-Hochberg adjusted
p-values, and the posterior null probabilities, respectively.
The second file contains a plot showing the values of \eqn{d_g} by gene No. on
the array.  The darker and larger the point, the smaller the adjusted p-value.
The gene with the smallest p-value is highlighted in green and its full name is
printed.
}
\usage{ lemmaPlots(outdir, mgq=0.99, titletext, modes) }

\arguments{
 \item{outdir}{The directory containing the RData file with the saved estimates
and data from a previous invocation of lemma(). Default=`OUT' } 
 \item{mgq}{The quantile used for plotting the \eqn{m_g} histogram.  Genes with 
\eqn{m_g} values outside the mgq quantile are excluded from the plot.
Default=0.99}
 \item{titletext}{The title of the plot}
 \item{modes}{The number of assumed (null + nonnull) groups. Can be either 2 or
3 (default=3). }
}

\details{
This function produces five diagnostics plots in two files (in postscript and
pdf formats, saved in the \code{outdir} directory).
The file LEMMAplots1 contains four plots:

(a) is the fitted distribution for the histogram of \eqn{d_g},

(b) is the fitted distribution for the histogram of \eqn{m_g},

(c) and (d) are `volcano plot' that show the distribution of the
Benjamini-Hochberg adjusted p-values, and the posterior null probabilities,
respectively.

The file LEMMAplots2 contains a plot showing the values of \eqn{d_g} by gene
number on the array.  The darker and larger the point, the smaller the
adjusted p-value.
The gene with the smallest p-value is highlighted in green and its full name is
printed.

Note - we are using ps2pdf to convert the postscript format to pdf. It is
assumed to be installed on the computer.
}

\keyword{ documentation }


