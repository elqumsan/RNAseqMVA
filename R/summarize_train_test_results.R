#' @title Summarize results of several training-testing experiments.
#' @author Jacques.van-Helden\@univ-amu.fr
#' @description Extract the results of a list containing training-testing
#' results produced by the function IterateTrainingTesting().
#' @param experimentList a list of results returned by IterateTrainingTesting()
#' @param verbose=1 level of verbosity
#'
#' @examples
#'
#' recountID <- "SRP042620"
#' summaryTable <- SummarizeTrainTestResults(
#'   experimentList = train.test.results.knn.k.values[[recountID]])
#'
#' @export
SummarizeTrainTestResults <- function(experimentList,
                                      experimentLabels = names(experimentList),
                                      verbose = 1) {
  if (!is(object = experimentList, class2 = "list")) {
    stop("SummarizeTrainTestResults() requires a list as experimentList argument")
  }

  ## Collect experiment IDs
  experiment.names <- names(experimentList)

  ## Check that there is at least one experiment
  if (length(experiment.names) < 1) {
    message("SummarizeTrainTestResults()\tempty list of results")
    return()
  }

  if (verbose >= 1) {
    message("\tSummarizing results from ", length(experiment.names), " training-testing experiments")
    if (verbose >= 2) {
      message("\tExperiment IDs\n\t\t", paste(collapse = "\n\t\t", experiment.names))
    }
  }


  ## Instantiate data frame to store testing and training error rates
  testing.error.rates <- data.frame()
  training.error.rates <- data.frame()

  ## Check that experiment labels are consistent with experiemnt list
  if (length(experiment.names) != length(experimentLabels)) {
    stop("length of experimentLabels (", length(experimentLabels),
         ") is inconsistent with length of experimentList (", length(experimentList),").")
  }

  ## Make sure the output directory exists for storing the plots
  if (is.null(experiment.names) ||
      length(names(experimentList)) < 1) {
    stop("Not a single experiment for all variable composition")
  }

  ## Collect all error rates in a data frame with 1 column per experiment and
  ## 1 row per iteration
  exp <- 1
  for (exp in 1:length(experiment.names)) {
    experimentName <- experiment.names[exp]
    exp.result <- experimentList[[experimentName]]
    if (exp == 1) {
      testing.error.rates <- data.frame(exp.result$testing.error.rate)
      training.error.rates <- data.frame(exp.result$training.error.rate)
    } else {
      testing.error.rates <-
        cbind(testing.error.rates, exp.result$testing.error.rate)
      training.error.rates <-
        cbind(training.error.rates, exp.result$training.error.rate)
    }
  } # end iterate the all variables experiment
  # dim(testing.error.rates)

  colnames(testing.error.rates) <- experimentLabels
  colnames(training.error.rates) <- experimentLabels

  #### Compute summary statistics per experiment ####
  testing.MER.summary <- data.frame(
    mean = apply(testing.error.rates, 2, mean),
    min = apply(testing.error.rates, 2, min),
    Q1 = apply(testing.error.rates, 2, quantile, prob = 0.25),
    median = apply(testing.error.rates, 2, median),
    Q2 = apply(testing.error.rates, 2, quantile, prob = 0.75),
    max = apply(testing.error.rates, 2, max)
  )
  # x <- data.frame(testing = testing.MER.summary$mean,
  #                 training = training.MER.summary$mean)
  # rownames(x) <- rownames(testing.MER.summary)
  # barplot(as.matrix(x), beside = TRUE)

  training.MER.summary <- data.frame(
    mean = apply(training.error.rates, 2, mean),
    min = apply(training.error.rates, 2, min),
    Q1 = apply(training.error.rates, 2, quantile, prob = 0.25),
    median = apply(training.error.rates, 2, median),
    Q2 = apply(training.error.rates, 2, quantile, prob = 0.75),
    max = apply(training.error.rates, 2, max)
  )


  #### Create result object ####
  result <- list(
    experiment.names = experiment.names,
    experimentLabels = experimentLabels,
    testing.error.rates = testing.error.rates,
    testing.MER.summary = testing.MER.summary,
    training.error.rates = training.error.rates,
    training.MER.summary = training.MER.summary)

  return(result)
}

