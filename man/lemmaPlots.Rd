\name{lemmaPlots}
\alias{lemmaPlots}
\docType{methods}
\title{LEMMA - diagnostics plots}
\description{
This program produces two diagnostics plots for the two fitting procedures in the 
LEMMA model.
}
\usage{ lemmaPlots(outdir, mgq=0.99, titletext, modes) }

\arguments{
 \item{outdir}{The directory containing the RData file with the saved estimates and
data from a previous invocation of lemma(). Default=`OUT' } 
 \item{mgq}{The quantile used for plotting the \eqn{m_g} histogram.  Genes with 
\eqn{m_g} values outside the mgq quantile are excluded from the plot.  Default=0.99}
 \item{titletext}{The title of the plot}
 \item{modes}{The number of assumed (null + nonnull) groups. Can be either 2 or 3
(default=3). }
}

\details{
The function lemmaPlots generates two plots on one page, and saves them in postscript and pdf formats under the \code{outdir} directory.

(a) a histogram of the \eqn{d_g} statistics and the fitted curves for the null and nonnull groups, as well as the mixture (using the average error variance).

(b) a histogram of \eqn{m_g}, and the fitted distribution, based on the estimates for alpha and beta.

Note - we are using ps2pdf to convert the postscript format to pdf. It is assumed to be installed on the computer.
}

\keyword{ documentation }


