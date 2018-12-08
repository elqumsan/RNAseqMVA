###############################################################################
#' @title Function for boxplot the error rate
#' @author Mustafa AbuElqumsan and Jacques van Helden
#' @description this function takes the experiment List from the error rate that are collected from the experiment.
#'
#' @param experimentList this is list of experiment where each cell of list error rate of the single experiment with special parameters
#' @param classifier such are represent which classifier have been used in the analysis
#' @param main is the main title of the boxplot
#' @param expMisclassificationRate user-specified value for the expected misclassification rate
#' @param expColor is colore each box plot that represent the error rate for each experiment with actual class lables.
#' @param permColor is the coler for the box plot that represent error rate for each experiment with permuted calss lables.
#' @param expLegend is the legend for the real class lable train/test experiment.
#' @param permLegend legend for the IterateTrainingTesting with permuted class lables.
#' @param legend.place location for the legend, which is passed to legend()
#'
#' @examples
#' compareExperiments(experimentList = train.test.results)
##################################################################################
## Gather error rates for each K Value in the experiments
# and draw a box plot
##################################################################################
#' @export
ErrorRateBoxPlot <- function(experimentList,
                             classifier,
                             expMisclassificationRate = NULL,
                             horizontal = TRUE,
                             fig.height = 10,
                             fig.width = 4 + 0.2*length(experimentList),
                             cex.axis = 0.8,
                             experimentLabels = names(experimentList),
                             main = paste(sep = "", parameters$recountID, "; ", classifier,
                                          "\n", parameters$iterations, "iterations"),
                             boxplotFile = NULL,
                             expColor = "#00BBFF",
                             permColor = "grey",
                             expLegend = "Train/test",
                             permLegend = "Permuted labels",
                             legend.place = "topright",
                             ....) {


  experimentSummary <- SummarizeTrainTestResults(experimentList = experimentList, experimentLabels = experimentLabels)

  testing.error.rates <- experimentSummary$testing.error.rates


  # ## Get experiment names to access the elements of experiment list
  # experiment.names <- names(experimentList)
  #
  # ## Check that experiment labels are consistent with experiemnt list
  # if (length(experiment.names) != length(experimentLabels)) {
  #   stop("length of experimentLabels (", length(experimentLabels),
  #        ") is inconsistent with length of experimentList (", length(experimentList),").")
  # }
  #
  # ## Make sure the output directory exists for storing the plots
  # if (is.null(experiment.names) ||
  #     length(names(experimentList)) < 1) {
  #   stop("Not a single experiment for all variable composition")
  # }
  #
  # ## Collect all error rates in a data frame with 1 column per experiment and
  # ## 1 row per iteration
  # exp <- 1
  # for (exp in 1:length(experiment.names)) {
  #   experimentName <- experiment.names[exp]
  #   exp.result <- experimentList[[experimentName]]
  #   if (exp == 1) {
  #     testing.error.rates <- data.frame(exp.result$testing.error.rate)
  #   } else {
  #     testing.error.rates <-
  #       cbind(testing.error.rates, exp.result$testing.error.rate)
  #   }
  # } # end iterate the all variables experiment
  # # dim(testing.error.rates)
  #
  # colnames(testing.error.rates) <- experimentLabels
  # rownames(testing.error.rates) <- 1:nrow(testing.error.rates)

  if (!is.null(boxplotFile)) {
    message.with.time("Boxplot file:\t", boxplotFile)
    if (horizontal) {
      pdf(file = boxplotFile, width = fig.height, height = fig.width)
    } else {
      pdf(file = boxplotFile, width = fig.width, height = fig.height)
    }
  }

  ## Define parameters for the boxplot
  save.margins <- par("mar")
  labelMargin <- (2 + max(nchar(experimentLabels)) * 0.5 * cex.axis)

  ## Define colors for experiments with actual data and permutation tests
  testTable.colors <- rep(x = expColor, length.out = ncol(testing.error.rates))
  names(testTable.colors) <- colnames(testing.error.rates)
  permTestExperiments <- grep(names(testTable.colors), pattern = "permLabels")
  if (length(permTestExperiments > 1)) {
    testTable.colors[permTestExperiments] <- permColor
  }

  ## Estimate the expected misclassification rate by computing the
  ## mean misclassification over all experiments with permuted labels
  ## Mean misclassification rate for all the label-permuted tests
  permLabelMisclassificationRate <- mean(apply(testing.error.rates[permTestExperiments], 1, mean))

  if (horizontal) {
    par(mar = c(5, labelMargin, 5, 1))
    ## Draw the box plot
    boxplot(testing.error.rates[, ncol(testing.error.rates):1],
            horizontal = TRUE,
            xlab = "Misclassification rate",
            ylim = c(0,1.3), ## Note: this actuallly corresponds to X limits with horizontal option
            # xlab = experimentLabels,
            main = main,
            las = 1 , cex.axis = cex.axis,
            col = rev(testTable.colors)
    )
    ## Draw horizontal grid
    abline(v = seq(from = 0, to = 1, by = 0.1), col = "grey", lty = "solid")

    ## Expected misclassification rate
    abline(v = permLabelMisclassificationRate, col = "blue", lwd = 1, lty = "dashed")
    if (!is.null(expMisclassificationRate)) {
      abline(v = expMisclassificationRate, col = "red", lwd = 3 , lty = "dotted")
    }


  } else {
    par(mar = c(labelMargin, 5, 5, 1))

    ## Draw the box plot
    boxplot(testing.error.rates,
            horizontal = FALSE,
            ylab = "Misclassification rate",
            ylim = c(0,1.1), ## Leave place for the legend
            # xlab = experimentLabels,
            main = main,
            las = 2 , cex.axis = cex.axis,
            col = testTable.colors
    )

    ## Draw horizontal grid
    abline(h = seq(from = 0, to = 1, by = 0.1), col = "grey", lty = "solid")
    #abline(h=seq(from = 0, to = 1, by= 0.05 ), col="", lty="solid")

    ## Draw a horizontal line with the misclassification
    ## rate that would be expected at random
    abline(h = permLabelMisclassificationRate, col = "blue", lwd = 1, lty = "dashed")
    if (!is.null(expMisclassificationRate)) {
      abline(h = expMisclassificationRate, col = "red", lwd = 3 , lty = "dotted")
    }
  }



  ## Draw a line with the theoretically computed random expected misclassification rate.
  ## In principle, this should fit the mean of the permutation tests.
  ## TO BE ONE BUT WE NEED TO EITHER PASS THE DATASET OR ASSOCIATE THIS METHOD TO THE CLASS.
  # abline(h = , col="red", lwd=3 , lty= 3)

  ## Plot legend
  legend(legend.place, lwd = 4,
         bty = "o", bg = "white",
         legend = c("Actual data", "Permuted Labels"),
         col = c(expColor, permColor),
         cex = 0.8, pch = 0.2)
  par(mar = save.margins)

  if (!is.null(boxplotFile)) {
    silence <- dev.off(); rm(silence)
  }

}

