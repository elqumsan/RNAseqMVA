% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tune_classifiers.R
\name{plotKNNtuning}
\alias{plotKNNtuning}
\title{draw plots to highlight the results of KNN tuning}
\usage{
plotKNNtuning(tuningResult, pdfFile = NULL, ...)
}
\arguments{
\item{tuningResult}{result of TuneClassifiers() run with option tuneKNN=TRUE, and where a vector of several values was provided for k}

\item{...}{additional parameters are passed to plot()}

\item{pdfFile=NULL}{if not null, the plot is exported in pdf format the specified file. Else it is printed out on the current device.}
}
\value{
no return field
}
\description{
draw a series of plots highlighting the impact of the number of neighbours (k) on the performances (miclassification error rate) of KNN
}
\author{
Jacques van Helden
}
