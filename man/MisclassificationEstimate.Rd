% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misclassification_estimate.R
\name{MisclassificationEstimate}
\alias{MisclassificationEstimate}
\title{RNA-Seq classifier evaluation to assess the performance of the classifier}
\usage{
MisclassificationEstimate(
  dataset,
  iteration,
  classifier,
  permute = FALSE,
  stratifiedPermutation = FALSE
)
}
\arguments{
\item{dataset}{an object of class DataTableWithTrainTestSets}

\item{iteration}{current iteration number (the MisclassificationEstimate function is typically called iteratively) and it called number of resampling.}

\item{classifier}{is a type of the classifier.}

\item{permute=FALSE}{if TRUE, permute the calss labels to measure misclassifciation rate without relevant training}
}
\value{
\itemize{
\item dataset = it is data table one row for the individual (sample) and one culomn for the feature (gene).
\item classifier = what is the tested classifier, KNN, RF, SVM and LDA.
\item iteration = nomber of resampling.
\item testPredictedClasses = these are the test classes which are predicted from classifier.
\item trainPredictedClasses = these are the train classes whci are predicated from classifier.
}
}
\description{
this script evaluates and assesses the performance of the RNA-Seq classifier by
random resampling (random partitioning) estimation of the misclassification error rate.
}
\author{
Mustafa ABUELQUMSAN and Jacques van Helden
}
