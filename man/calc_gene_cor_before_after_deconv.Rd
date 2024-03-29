% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TME-deconvolution.R
\name{calc_gene_cor_before_after_deconv}
\alias{calc_gene_cor_before_after_deconv}
\title{calculate expression of loci to specific genes before and after TME deconvolution}
\usage{
calc_gene_cor_before_after_deconv(deconv_list, raw_meth, genes)
}
\arguments{
\item{deconv_list}{output of \code{deconv_TME}.}

\item{raw_meth}{intervals set with raw methylation}

\item{genes}{vector with gene names}
}
\value{
data frame with intervals set that are in raw_meth and normalized methylation and their correlations to the gene expression in the raw and normalized methylation datasets.
}
\description{
calculate expression of loci to specific genes before and after TME deconvolution
}
