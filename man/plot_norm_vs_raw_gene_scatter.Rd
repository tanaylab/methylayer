% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize-em-cross.R
\name{plot_norm_vs_raw_gene_scatter}
\alias{plot_norm_vs_raw_gene_scatter}
\title{Plot scatter of raw expression-methylation correaltion vs normalized for specific genes}
\usage{
plot_norm_vs_raw_gene_scatter(em_list_raw, em_list_norm, gene, color = "black")
}
\arguments{
\item{em_list_raw}{output of \code{em_cross_cor} on raw data}

\item{em_list_norm}{output of \code{em_cross_cor} on normalized data}

\item{gene}{name of the gene to plot}

\item{color}{points color}
}
\value{
ggplot object
}
\description{
Plot scatter of raw expression-methylation correaltion vs normalized for specific genes
}
