\name{printTopGenes}
\alias{printTopGenes}
\docType{methods}
\title{LEMMA - printing the most likely nonnull genes}
\description{
This program prints `topgenes' genes sorted in ascending order by their 
posterior null probability, or by their adjusted p-value. Use a positive
integer, or `all' to print all the genes, or `nonnull' to print only the genes
that are declared nonnull with the given local fdr or FDR cut-off values. 
It is executed by \code{lemma()}, but can be executed later separately, if the
user wants to print a different number of `topgenes', without running the EM
algorithm again.
}
\usage{
  printTopGenes(type, outdir, data0, data1, data2, geneid, genename,
     topgenes, titletext, cutoff, modes)
}
\arguments{
 \item{type}{Either `RR' or `FDR'.}
 \item{outdir}{The directory with the AllData.RData file from a previous
invocation of \code{lemma()}. The output files will be stored in the outdir
directory (overwriting the previous copy).}
 \item{data0}{Enter RRfdr0 if type=`RR', or pBH0 if type=`FDR'}
 \item{data1}{Enter RRfdr1 if type=`RR', or (dg-tau) if type=`FDR'}
 \item{data2}{Enter (dg-tau) if type=`RR', or c() if type=`FDR'}
 \item{geneid}{Use the geneid vector}
 \item{genename}{Use the genename vector}
 \item{topgenes}{the number of genes to print (sorted in ascending order by
local-fdr or adjusted p-value). Use `all' to print all the genes, or `nonnull'
to print only the genes that are declared nonnull with the given local fdr or
FDR cut-off values. Default=`nonnull'.}
 \item{titletext}{The title to be used in the output files.}
 \item{cutoff}{Enter locfdrcutoff if type=`RR', or fdrcutoff if type=`FDR'}
 \item{modes}{the number of assumed (null + nonnull) groups. Can be either
 2 or 3 (default=3)}
}

\keyword{ documentation }

