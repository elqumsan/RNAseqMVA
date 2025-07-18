% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/study_case.R
\name{RunedgeR.StudyCase}
\alias{RunedgeR.StudyCase}
\title{run edgeR on an object of class StudyCase}
\usage{
\method{RunedgeR}{StudyCase}(self)
}
\arguments{
\item{self}{object belong to StudyCase class}
}
\value{
a clone of the input StudyCase object with an added
DataTableWithTrainTestSets containing the log2-normalized counts data table
where features have been re-ordered by increasing adjusted p-value.
}
\description{
run edgeR on an object of class StudyCase to test differential expression on each feature of a data table, and order variables by increasing adjusted p-value.
}
\author{
Mustafa AbuElQumsan and Jacques van Helden
}
