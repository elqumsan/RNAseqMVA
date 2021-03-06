% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tune_classifiers.R
\name{TuneClassifiers}
\alias{TuneClassifiers}
\title{Tune the parameters of the different classifiers for a given study case}
\usage{
TuneClassifiers(
  studyCase,
  dataType = "TMM_log2_PC",
  tuneSVM = TRUE,
  svmCost = 1,
  svmDegree = 1:4,
  svmGamma = 4^(-2:4),
  svmCoef0 = 4^(-2:4),
  tuneRandomForest = TRUE,
  rfNtree = c(100, 200, 300, 500),
  rfMtry = c(50, 100, 200, 300, 500),
  tuneKNN = TRUE,
  knnK = c(1, 2, 3, 4, 5, 7, 10, 15),
  knnL = 0,
  plotResults = FALSE
)
}
\arguments{
\item{studyCase}{a studyCase object}

\item{dataType="TMM_log2_PC"}{data type to use for the tuning. Must be one of the data types included in the 'datasetsForTest' attribute of the studyCase object}

\item{tuneSVM=TRUE}{if TRUE, tune parameters for SVM}

\item{svmCost=1}{cost parameter values for SVM}

\item{svmDegree=1:4}{degree values for SVM polynomial kernel}

\item{svmGamma=4^(-2:4)}{gamma values for SVM kernels (ignored for linear kernel)}

\item{tuneRandomForest=TRUE}{if TRUE, tune parameters for Random Forest}

\item{rfNtree=c(100, 200, 300, 500)}{values to test for RF numbers of trees}

\item{rfMtry=c(50, 100, 200, 300, 500)}{values to test for RF mtry parameter}

\item{tuneKNN=TRUE}{if TRUE, tune parameters for KNN}

\item{knnK=c(1, 2, 3, 4, 5, 7, 10, 15)}{number of neighbours for KNN}

\item{knnL=0}{minimum vote for definite decision for KNN}

\item{plotResults=FALSE}{if TRUE, run plot()  on the tuned objects}
}
\value{
a list of parameters cloned from the studyCase, added with the optimal parameters identified here
}
\description{
Tune the parameters of the different classifiers for a given study case
}
\author{
Jacques van Helden
}
