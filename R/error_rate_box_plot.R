###############################################################################
#' @title Function for boxplot the error rate
#' @author Mustafa AbuElqumsan and Jacques van Helden
#' @description this function takes the experiment List from the error rate that are collected from the experiment.
#'
#' @param experimentList this is list of experiment where each cell of list error rate of the single experiment with special parameters
#' @param classifier name of the classifier (e.g. svm, knn, rf), used in the default main title
#' @param main main title for the boxplot
#' @param expMisclassificationRate user-specified value for the expected misclassification rate
#' @param medianLine=FALSE if TRUE, draw a line joining the medians
#' @param expColor is colore each box plot that represent the error rate for each experiment with actual class lables.
#' @param permColor is the coler for the box plot that represent error rate for each experiment with permuted calss lables.
#' @param expLegend is the legend for the real class label train/test experiment.
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
                             fig.width = 4 + 0.05 * length(experimentList),
                             cex.axis = 0.8,
                             experimentLabels = names(experimentList),
                             main = paste0(
                               parameters$short_label,
                               " (", parameters$recountID, ") ", parameters$feature, "s",
                               "\n", classifier, "; ", parameters$iterations, " iterations"),
                             boxplotFile = NULL,
                             medianLine = FALSE,
                             expColor = "#00BBFF",
                             permColor = "grey",
                             expLegend = "Train/test",
                             permLegend = "Permuted labels",
                             legend.place = "topright",
                             ....) {


  experimentSummary <- SummarizeTrainTestResults(experimentList = experimentList, experimentLabels = experimentLabels)

  testing.error.rates <- experimentSummary$testing.error.rates

  if (!is.null(boxplotFile)) {
    message.with.time("Boxplot file:\t", boxplotFile)
    if (horizontal) {
      pdf(file = boxplotFile, width = fig.height, height = fig.width)
    } else {
      pdf(file = boxplotFile, width = fig.width, height = fig.height)
    }
  }


  ## Define colors for experiments with actual data and permutation tests
  testTable.colors <- rep(x = expColor, length.out = ncol(testing.error.rates))
  names(testTable.colors) <- colnames(testing.error.rates)
  permTestExperiments <- grep(names(testTable.colors), pattern = "permLabels")
  dataColumns <- ncol(testing.error.rates):1

  ## Define parameters for the boxplot
  save.margins <- par("mar")
  labelMargin <- (2 + max(nchar(experimentLabels[dataColumns])) * 0.46 * cex.axis)

  if (length(permTestExperiments > 1)) {
    testTable.colors[permTestExperiments] <- permColor
    permColumns <- rev(permTestExperiments)
    dataColumns <- setdiff(dataColumns, permColumns)
    ## Estimate the expected misclassification rate by computing the
    ## mean misclassification over all experiments with permuted labels
    ## Mean misclassification rate for all the label-permuted tests
    permLabelMisclassificationRate <- mean(apply(testing.error.rates[permTestExperiments], 1, mean))
  }



  if (horizontal) {
    par(mar = c(5, labelMargin, 5, 1))

    ## Draw the box plot
    boxplot(testing.error.rates[, dataColumns],
            horizontal = TRUE,
            xlab = "Misclassification error rate (MER)",
            ylim = c(0,1.3), ## Note: this actuallly corresponds to X limits with horizontal option
            # xlab = e`xperimentLabels,
            main = main,
            las = 1 , cex.axis = cex.axis,
            col = testTable.colors[dataColumns]
    )

    ## Permuted labels
    if (length(permTestExperiments > 1)) {
      abline(v = permLabelMisclassificationRate, col = "blue", lwd = 1, lty = "dashed")
      boxplot(testing.error.rates[, permColumns],
              horizontal = TRUE,
              col = testTable.colors[permColumns],
              las = 1,
              names = NA,
              cex.axis = 0.01,
              add = TRUE
      )

      ## Draw horizontal grid
      abline(v = seq(from = 0, to = 1, by = 0.05), col = "grey", lty = "dotted")
      abline(v = seq(from = 0, to = 1, by = 0.1), col = "grey", lty = "solid")

    }
    if (!is.null(expMisclassificationRate)) {
      abline(v = expMisclassificationRate, col = "red", lwd = 3 , lty = "dotted")
    }

    ## Draw lines between medians
    if (medianLine) {
      lines(rev(experimentSummary$testing.MER.summary$median[dataColumns]), dataColumns, col = testTable.colors[dataColumns], lwd = 2)
      lines(rev(experimentSummary$testing.MER.summary$median[permColumns]), dataColumns, col = testTable.colors[permColumns], lwd = 2)
    }

  } else {
    par(mar = c(labelMargin, 5, 5, 1))

    ## Draw the box plot
    boxplot(testing.error.rates[, dataColumns],
            horizontal = FALSE,
            ylab = "Misclassification error rate (MER)",
            ylim = c(0,1.1), ## Leave place for the legend
            # xlab = experimentLabels,
            main = main,
            las = 2 , cex.axis = cex.axis,
            col = testTable.colors[dataColumns]
    )

    ## Boxplots of permuted labels
    if (length(permTestExperiments > 1)) {
      boxplot(testing.error.rates[, permColumns],
              horizontal = FALSE,
              col = testTable.colors[permColumns],
              las = 1, names = NA, cex.axis = 0.01,
              add = TRUE
      )
    }

    ## Draw grid
    abline(h = seq(from = 0, to = 1, by = 0.05), col = "grey", lty = "dotted")
    abline(h = seq(from = 0, to = 1, by = 0.1), col = "grey", lty = "solid")

    ## Draw lines between medians
    if (medianLine) {
      lines(dataColumns, rev(experimentSummary$testing.MER.summary$median[dataColumns]), col = testTable.colors[dataColumns], lwd = 2)
      lines(dataColumns, rev(experimentSummary$testing.MER.summary$median[permColumns]), col = testTable.colors[permColumns], lwd = 2)
    }


    ## Draw a horizontal line with the misclassification
    ## rate that would be expected at random
    if (length(permTestExperiments > 1)) {
      abline(h = permLabelMisclassificationRate, col = "blue", lwd = 1, lty = "dashed")
    }
    if (!is.null(expMisclassificationRate)) {
      abline(h = expMisclassificationRate, col = "red", lwd = 3 , lty = "dotted")
    }
  }



  ## Draw a line with the theoretically computed random expected misclassification rate.
  ## In principle, this should fit the mean of the permutation tests.
  ## TO BE ONE BUT WE NEED TO EITHER PASS THE DATASET OR ASSOCIATE THIS METHOD TO THE CLASS.
  # abline(h = , col="red", lwd=3 , lty= 3)

  ## Plot legend
  legend.text <- "Actual data"
  legend.col <- expColor
  legend.lty <- "solid"
  legend.lwd <- 4
  legend.pch <- 1
  if (length(permTestExperiments > 1)) {
    legend.text <- append(legend.text, c("Permuted Labels", "permuted MER"))
    legend.col <- append(legend.col, c(permColor, "blue"))
    legend.lty <- append(legend.lty, c("solid", "dotted"))
    legend.lwd <- append(legend.lwd, c(4, 2))
    legend.pch <- append(legend.pch, c(1, NA))
  }
  if (!is.null(expMisclassificationRate)) {
    legend.text <- append(legend.text, "expected MER")
    legend.col <- append(legend.col, "red")
    legend.lty <- append(legend.lty, "dotted")
    legend.lwd <- append(legend.lwd, 2)
    legend.pch <- append(legend.pch, NA)
  }

  legend(legend.place,
         bty = "o",
         bg = "white",
         legend = legend.text,
         col = legend.col,
         lty = legend.lty,
         lwd = legend.lwd,
         cex = 0.7, pch = legend.pch)
  par(mar = save.margins)

  if (!is.null(boxplotFile)) {
    silence <- dev.off(); rm(silence)
  }

}

