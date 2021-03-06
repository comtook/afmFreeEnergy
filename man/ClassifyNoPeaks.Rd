% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClassifyNoPeaks.R
\name{ClassifyNoPeaks}
\alias{ClassifyNoPeaks}
\title{The main function classifies data into five classes
       according to number of binding points}
\usage{
ClassifyNoPeaks(datDir, desDir, Alpha = 0.5, LOD = 3, nlsEr = 4,
  mulBase = 2)
}
\arguments{
\item{datDir}{data directory}

\item{desDir}{result directory}

\item{Alpha}{percent of data points in bulk region}

\item{LOD}{user given level of noise in standard deviation unit}

\item{nlsEr}{RMSD tolerence for nls fit in standard deviation unit}

\item{mulBase}{user given level of baseline perturbation in standard deviation unit}
}
\description{
This function allows you to classify F-Z curves into five classes.
}
