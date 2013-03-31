\name{lemma}
\alias{lemma}
\docType{methods}
\title{LEMMA - the main program}
\description{
This program implements the EM algorithm for the LEMMA model.
LEMMA fits a linear mixed-model for normalized
microarray data analysis. See the reference to the LEMMA paper, which contains
the underlying model and the theory.  For a short summary of the model, see

http://www.stat.cornell.edu/lemma/docs/LEMMAsummary.pdf.

This version supports two treatment groups
and either a two-way classification (null and nonnull genes),
or a three-way classification: null genes, for which statistically
there is no difference in expression between the two treatment groups;
nonnull group 1 - genes that are significantly more expressed in treatment
group 1 than in treatment group 2;
and nonnull group 2 - genes that are significantly more expressed in treatment
group 2 than in treatment group 1.
This version also supports within-group analysis (for example, in paired
experiments).

The program runs on both Windows and Linux.
}

\details{
The input should consist of a data frame with G rows (one row per gene),
and have the following structure:
\itemize{
\item{geneid}{ - the first column must contain a vector of G gene IDs}
\item{genename}{ - the second column must contain a vector of G gene names}
\item{Y1}{ - columns 3,...,(2+n1) contain normalized data for n1 replicates of
subjects from treatment group 1.  The columns must be named Y1_m where
m=1,...,n1.}
\item{Y2}{ - columns (n2+2),...,(2+n1+n2) contain normalized data for n2
replicates of subjects from treatment group 2.  The columns must be named
Y2_m where m=1,...,n2. Note that these columns are only used in
between-treatment analysis.  In the case of within-group analysis, these
columns should not be used. }
}

In this version n1 and n2 do not have to be the same, but all the rows in Y1
have to have n1 elements, and all the rows in Y2 have to have n2 elements.

All of the parameter estimates, plots, and gene lists will be saved under the
user-specified outdir directory. In particular, this directory will contain
the following files:
\itemize{
\item{\code{log.txt}}{ - reporting the total number of genes, sample sizes,
\eqn{mean(d_g), sd(d_g), mean(m_g),} \eqn{sd(m_g)}, estimates of the shape and scale
parameters of the assumed inverse gamma prior for the error variance. It also
contains the mean and variance of the fitted error variance distribution
(they should be close to the sample mean and variance based on the observed
\eqn{m_g}).
Estimates for \eqn{\tau, \psi, \sigma^2_\psi, p_1} and \eqn{p_2} are also
included in this log file, as well as the number of nonnulls genes detected
using the user-provided local fdr cutoff, and the FDR threshold. Note that in
this version parameter estimates are provided with their standard errors
(based on the Fisher Information matrix for the error variance parameter,
and the Oakes method for the mean and mixture parameters).
Any convergence problems in the EM algorithm are reported in this file.}
\item{\code{resultsRR.txt}}{ - contains a list of genes sorted by their
posterior null probability. This file also contains the estimated posterior
probabilities for a gene being more expressed in treatment group 1 than in
treatment group 2 (and vice versa). It also contains the gene effect
(\eqn{d_g-\tau}). }
\item{\code{resultsFDR.txt}}{ - contains a list of genes sorted by their
BH-adjusted p-values. The file also contains the gene effect (\eqn{d_g-\tau}),
and the sign of the gene effect which can be used to determine if a (nonnull)
gene is more expressed in treatment group 1 than in treatment group 2 (or vice
versa).}
\item{\code{results.csv}}{ - contains a list of all the genes, their
BH-adjusted p-values, the test statistic (\eqn{d_g}), and the posterior
probabilities.}
\item{\code{AllData.RData}}{ - contains the following elements: dg, mg, n1, n2,
f, G, RRfdr0, RRfdr1, RRfdr2, alpha_hat, beta_hat, sig2eb, tau, psi, sig2psi,
p0, p1, p2, pBH0}}

Note:\itemize{
\item{f}{ = \eqn{n_1+n_2-2} (the degrees of freedom for the mean square error,
\eqn{m_g}) }
\item{sig2eb}{ = a vector of length G containing the posterior mode of the
error variances, given \eqn{m_g} }
\item{pBH0}{ = a vector of length G containing the BH-adjusted p-value of genes}
\item{RRfdr0}{ = a vector of length G containing the posterior probabilities
of genes being null}
\item{RRfdr1}{ = a vector of length G containing the posterior probabilities
of genes being nonnull and more expressed in treatment group 1 than in
treatment group 2}
\item{RRfdr2}{ = a vector of length G containing the posterior probabilities
of genes being nonnull and more expressed in treatment group 2 than in
treatment group 1}
}
}

\usage{ lemma(dataframe, locfdrcutoff=0.2, fdrcutoff=0.2, mgq=1,
    titletext="", outdir=tempdir(), topgenes="nonnull", tol=1e-6,
    maxIts=50000, modes=3, plots=TRUE, saveascsv=TRUE,
    ErrVarEst="MLE") }

\arguments{
 \item{dataframe}{The data frame file containing the normalized data (the
first column contains the gene IDs, the second column contains the gene names,
the next n1 columns contain the normalized data for treatment group 1 (Y1),
and in the two-treatment case, the last n2 columns 
contain the normalized data for treatment group 2 (Y2).}
 \item{locfdrcutoff}{the local fdr cutoff value for detecting nonnull genes.
Default=0.2.}
 \item{fdrcutoff}{the level for the FDR procedure for detecting nonnull genes.
Default=0.2.}
 \item{mgq}{The quantile used for eliminating genes with extreme \eqn{m_g}
values.  Default=1.}
 \item{titletext}{The title to be used in the output files.}
 \item{outdir}{the output directory (use / as a directory separator, even on
Windows).  Default=tempdir().}
 \item{topgenes}{the number of genes to print (sorted in ascending order by
local-fdr or adjusted p-value). Use `all' to print all the genes, or `nonnull'
to print only the genes that are declared nonnull with the given local fdr or
FDR cut-off values. Default=`nonnull'.}
 \item{tol}{The tolerance level for determining the convergence of the EM
algorithm.\cr Default=1e-6.}
 \item{maxIts}{The maximum number of EM iterations.  Default=50000.}
 \item{modes}{the number of assumed (null + nonnull) groups. Can be either 2
or 3 (default=3). }
 \item{plots}{logical - if TRUE, produce the diagnostics plots in ps and pdf
formats (default=TRUE). }
 \item{saveascsv}{logical - if TRUE, produce the comma-separated file with
all the genes and their adjusted p-values and posterior probabilities
(default=TRUE). }
 \item{ErrVarEst}{The estimation mathod for the error variance paremeters.
Can be either MLE for maximum likelihood estimation, or MM
for Method of Moments estimation (default=MLE).}
}


\examples{
\dontrun{
lemma(apoai,titletext="APO-AI, Callow et al (2000)", plots=F)
lemmaPlots(tempdir(),mgq=0.99, titletext="APO-AI (Callow et al., 2000)")

lemma(simdata,titletext="Simulated data")

# Similarly, if the user wants to use the 2-way classification:
lemma(apoai,titletext="APO-AI, Callow et al (2000)", modes=2, plots=F)
lemmaPlots(tempdir(),mgq=0.99,titletext="APO-AI (Callow et al., 2000)",
      ErrVarEst="MM",modes=2)
}
}
  
\keyword{ documentation }

