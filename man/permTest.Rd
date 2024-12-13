% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoSig.R
\name{permTest}
\alias{permTest}
\title{Permutation on samples to extract the genes}
\usage{
permTest(
  t_coexp,
  meta,
  varFeatures,
  featureOfInterest,
  featureOfCategory,
  prefix,
  nsamples_per_group = 5
)
}
\arguments{
\item{meta}{metadata for samples}

\item{varFeatures}{gene names}

\item{featureOfInterest}{column name of the feature of interest in the metadata assuming the features are levels or numerical values}

\item{nsamples}{number of samples}

\item{cor_method}{correlation method}

\item{zthres}{threshold for z-score of loading values to identify significant driving genes per TC}
}
\description{
Permutation on samples to extract the genes
}
