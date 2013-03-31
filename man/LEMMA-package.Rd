\name{LEMMA}
\alias{LEMMA}
\docType{package}
\title{Laplace approximated EM Microarray Analysis}
\description{
LEMMA fits a linear mixed model to normalized microarray data.
See the complete LEMMA paper (on the lemma web site) which contains
the underlying model and the theory.

http://www.stat.cornell.edu/lemma/docs/LEMMAsummary.pdf contains a short
summary of the \cr model.

This version supports either two treatment groups or within-group analysis (for
example, in paired tests).  It provides either a two-way classification
(null and nonnull genes), or a three-way classification: null genes, for which
statistically there is no difference in expression between the two treatment
groups; nonnull group 1 - genes that are significantly more expressed in
treatment group 1 than in treatment group 2;
and nonnull group 2 - genes that are significantly more expressed in treatment
group 2 than in treatment group 1.

The program runs on both Windows and Linux.
}
\details{
The input should consist of a data frame with G rows, and have the following
structure:
\itemize{
\item{geneid}{ - the first column must contain a vector of G gene IDs}
\item{genename}{ - the second column must contain a vector of G gene names}
\item{Y1}{ - columns 3,...,(2+n1) contain normalized data for n1 replicates of
subjects from treatment group 1.  The columns must be named Y1_m where
m=1,...,n1.}
\item{Y2}{ - columns (n2+2),...,(2+n1+n2) contain normalized data for n2
replicates of subjects from treatment group 2.  The columns must be named
Y2_m where m=1,...,n2.  Note that these columns should only be used in the
two-groups case.  In the case of within-group analysis, only Y1_m columns
should be in the data frame (as well as the gene id and gene name columns.)} 
}

In this version n1 and n2 do not have to be the same, but all the rows in Y1
have to have n1 elements, and all the rows in Y2 have to have n2 elements.
The program also uses the following variables when the user invokes the
\code{lemma} function:
\code{outdir}, \code{locfdrcutoff}, \code{fdrcutoff}, \code{topgenes},
\code{titletext}, \code{mgq}, \code{tol}, \code{maxIts}, \code{modes},
\code{plots}, \code{saveascsv}, \code{ErrVarEst}.

All of the parameter estimates, plots, and gene lists will be saved under the
user-specified outdir directory. In particular, this directory will contain
the following files:
\itemize{
\item{\code{log.txt}}{ - reporting the total number of genes, sample sizes,
 \eqn{mean(d_g), sd(d_g), mean(m_g),} \eqn{sd(m_g)}, estimates of the shape
 and scale parameters of the assumed inverse gamma prior for the error variance.
 It also contains the mean and variance of the fitted error variance
 distribution (they should be close to the sample mean and variance based on
 the observed \eqn{m_g}).
 Estimates for \eqn{\tau, \psi, \sigma^2_\psi, p_1} and \eqn{p_2}, as well
 as standard errors for all the parameters are also included in this log file.
 Also logged are the number of nonnulls genes detected using the user-provided
 local fdr and the FDR thresholds. Any convergence problems
 in the EM algorithm are reported in this file.}
\item{\code{resultsRR.txt}}{ - contains a list of genes sorted by their
posterior null probability. This file also contains the estimated posterior
probabilities for a gene being more expressed in treatment group 1 than in
treatment group 2 (and vice versa). It also contains the gene effect
(\eqn{d_g-\tau} ).}
\item{\code{resultsFDR.txt}}{ - contains a list of genes sorted by their
BH-adjusted p-values. The file also contains the gene effect (\eqn{d_g-\tau} ),
and the sign of the gene effect which can be used to determine if a (nonnull)
gene is more expressed in treatment group 1 than in treatment group 2 (or vice versa).}
\item{\code{results.csv}}{ - contains a list of all the genes, their
BH-adjusted p-values, the test statistic (\eqn{d_g}), and the posterior
probabilities.}
\item{\code{AllData.RData}}{ - contains the following elements:
dg, mg, n1, n2, f, G, RRfdr0, RRfdr1, RRfdr2, alpha_hat, beta_hat, sig2eb,
tau, psi, sig2psi, p0, p1, p2, pBH0}
}

Note:
\itemize{
\item{f}{ = \eqn{n_1+n_2-2} (the degrees of freedom for the mean square error, \eqn{m_g} )}
\item{sig2eb}{ = a vector of length G containing the posterior mode of the error variances,
 given \eqn{m_g} }
\item{pBH0}{ = a vector of length G containing the BH-adjusted p-value of genes}
\item{RRfdr0}{ = a vector of length G containing the posterior probabilities of genes being null}
\item{RRfdr1}{ = a vector of length G containing the posterior probabilities of genes being nonnull and more
expressed in treatment group 1 than in treatment group 2}
\item{RRfdr2}{ = a vector of length G containing the posterior probabilities of genes being
 nonnull and more expressed in treatment group 2 than in treatment group 1}
}
}
\author{
 Bar, H.Y. \email{hyb2@cornell.edu}, Schifano, E.D. \email{eds27@cornell.edu}
}
\references{
 Bar, H.Y., Booth, J.G., Schifano, E.D., Wells, M.T., (2010).
 Laplace approximated EM Microarray Analysis: an empirical Bayes approach for
 comparative microarray experiments.

 http://www.stat.cornell.edu/lemma/docs/lemma.pdf
}
\keyword{ documentation }
\seealso{
 Read \link{lemma} to see how to execute the program.

 Use \link{lemmaPlots} to produce diagnostics plots.

 Use \link{printTopGenes} to produce a list of genes sorted by their adjusted
p-values or by their posterior null probabilities.
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
\dontrun{
lemma(apoai,titletext="APO-AI, Callow et al (2000)", plots=F)
lemmaPlots(tempdir(),mgq=0.99, titletext="APO-AI (Callow et al., 2000)")

lemma(simdata,titletext="Simulated data")

# Similarly, if the user wants to use the 2-way classification:
lemma(apoai,titletext="APO-AI, Callow et al (2000)",outdir=tempdir(),
  modes=2, plots=F)
lemmaPlots(tempdir(),mgq=0.99,titletext="APO-AI (Callow et al., 2000)",
      modes=2)
}
}
