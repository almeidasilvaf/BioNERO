% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gcn_inference.R
\name{get_hubs_gcn}
\alias{get_hubs_gcn}
\title{Get GCN hubs}
\usage{
get_hubs_gcn(exp, net)
}
\arguments{
\item{exp}{A gene expression data frame with genes in row names and
samples in column names or a `SummarizedExperiment` object.}

\item{net}{List object returned by \code{exp2gcn}.}
}
\value{
Data frame containing gene IDs, modules and
intramodular connectivity of all hubs.
}
\description{
Get GCN hubs
}
\examples{
data(filt.se)
gcn <- exp2gcn(filt.se, SFTpower = 18, cor_method = "pearson")
hubs <- get_hubs_gcn(filt.se, gcn)
}
\seealso{
\code{\link[WGCNA]{signedKME}}
}
\author{
Fabricio Almeida-Silva
}
