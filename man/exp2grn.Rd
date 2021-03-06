% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grn_inference.R
\name{exp2grn}
\alias{exp2grn}
\title{Infer gene regulatory network from expression data}
\usage{
exp2grn(
  exp,
  regulators = NULL,
  eps = 0,
  estimator_aracne = "spearman",
  estimator_clr = "pearson",
  remove_zero = TRUE,
  nsplit = 10,
  ...
)
}
\arguments{
\item{exp}{A gene expression data frame with genes in row names and
samples in column names or a `SummarizedExperiment` object.}

\item{regulators}{A character vector of regulators
(e.g., transcription factors or miRNAs). All regulators must be
included in `exp`.}

\item{eps}{Numeric value indicating the threshold used when
removing an edge: for each triplet of nodes (i,j,k), the weakest edge,
say (ij), is removed if its weight is below min{(ik),(jk)} - eps. Default: 0.}

\item{estimator_aracne}{Entropy estimator to be used in ARACNE inference.
One of "mi.empirical", "mi.mm", "mi.shrink", "mi.sg", "pearson", "spearman",
or "kendall". Default: "spearman".}

\item{estimator_clr}{Entropy estimator to be used in CLR inference.
One of "mi.empirical", "mi.mm", "mi.shrink", "mi.sg", "pearson",
"spearman", or "kendall". Default: "pearson".}

\item{remove_zero}{Logical indicating whether to remove edges whose
weight is exactly zero. Zero values indicate edges that were
removed by ARACNE. Default: TRUE.}

\item{nsplit}{Number of groups in which the edge list will be split.
Default: 10.}

\item{...}{Additional arguments passed to `GENIE3::GENIE3()`.}
}
\value{
A filtered edge list with regulators in the first column and
targets in the second column.
}
\description{
Infer gene regulatory network from expression data
}
\details{
This function infers GRNs with ARACNE, GENIE3 and CLR, ranks correlation
weights for each GRN and calculates the average rank for each edge.
Then, the resulting GRN is filtered to keep the top n edges that lead
to the optimal scale-free topology fit.
}
\examples{
data(filt.se)
tfs <- sample(rownames(filt.se), size=50, replace=FALSE)
# Test with small number of trees for demonstration purpose
grn <- exp2grn(filt.se, regulators = tfs, nTrees=2, nsplit=2)
}
