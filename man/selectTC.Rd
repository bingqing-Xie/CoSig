% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoSig.R
\name{selectTC}
\alias{selectTC}
\title{Extract the tensor component (TC) that correlated with the feature of interest}
\usage{
selectTC(
  meta,
  featureOfInterest,
  nsamples,
  cp_decomp,
  varFeatures,
  cor_method = "spearman"
)
}
\arguments{
\item{meta}{metadata for samples}

\item{featureOfInterest}{column name of the feature of interest in the metadata: assuming the features are levels or numerical values}

\item{nsamples}{number of samples}

\item{cp_decomp}{decomposed matrices from HOSVD}

\item{varFeatures}{gene names}

\item{cor_method}{correlation method}
}
\value{
A list of selected TCs
}
\description{
Extract the tensor component (TC) that correlated with the feature of interest
}
