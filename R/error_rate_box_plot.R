###############################################################################
#' @title Function for boxplot the error rate
#' @author Mustafa AbuElqumsan and Jacques van Helden
#' @description this function takes the experiment List from the error rate that are collected from the experiment.
#'
#' @param experimentList this is list of experiment where each cell of list error rate of the single experiment with special parameters
#' @param classifier such are represent which classifier have been used in the analysis
#' @param main is the main title of the boxplot
#' @param expColor is colore each box plot that represent the error rate for each experiment with actual class lables.
#' @param permColor is the coler for the box plot that represent error rate for each experiment with permuted calss lables.
#' @param expLegend is the legend for the real class lable train/test experiment.
#' @param  permLegend legend for the IterateTrainingTesting with permuted class lables.
#' @param dataType data type(s), used to build the file prefix
#' @param variablesType variable types (e.g. all, DEG, ...), used to build the file prefix
#'
#' @examples
#' compareExperiments(experimentList = train.test.results)
##################################################################################
## Gather error rates for each K Value in the experiments
# and draw a box plot
##################################################################################
#' @export
ErrorRateBoxPlot <- function (
  experimentList,
  classifier,
  data.type,
  fig.height = 10,
  fig.width = 4 + 0.2*length(experimentList),
  experimentLabels = names(experimentList),
  main = paste(sep = "", parameters$recountID, "; ", classifier,
               "\n", data.type,
               #"; ", variablesType,
               "; ", parameters$iterations, "iterations"),
  boxplotFile = NULL,
  expColor = "#00BBFF",
  permColor = "grey",
  expLegend = "Train/test",
  permLegend = "Permuted labels",
  ....) {


  testing.error.rates <- data.frame()

  ## Get experiment names to access the elements of experiment list
  experiment.names <- names(experimentList)

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


  exp <-1
  for (exp in 1:length(experiment.names)) {
    experimentName <- experiment.names[exp]
    exp.result <- experimentList[[experimentName]]
    if (exp == 1) {
      testing.error.rates <- data.frame(exp.result$testing.error.rate)
    } else {
      testing.error.rates <-
        cbind(testing.error.rates,
              exp.result$testing.error.rate)
    }
  } # end iterate the all variables experiment

  # dim(testing.error.rates)
  colnames(testing.error.rates) <- experimentLabels
  rownames(testing.error.rates) <- 1:nrow(testing.error.rates)

  if (!is.null(boxplotFile)) {
    message.with.time("Boxplot file:\t", boxplotFile)
    pdf(file = boxplotFile, width = fig.width, height = fig.height)
  }

  ## Define parameters for the boxplot
  save.margins <- par("mar")
  labelMargin <- (2 + max(nchar(experimentLabels)) * 0.45)
  par(mar=c(labelMargin, 5, 5, 1))

  ## Define colors for experiments with actual data and permutation tests
  testTable.colors <- rep(x = expColor, length.out=ncol(testing.error.rates))
  names(testTable.colors) <- colnames(testing.error.rates)
  permTestExperiments <- grep(names(testTable.colors), pattern = "permLabels")
  if (length(permTestExperiments > 1)) {
    testTable.colors[permTestExperiments] <- permColor
  }


  ## Draw the box plot
  boxplot(testing.error.rates,
          horizontal = FALSE ,
          ylab = "Misclassification rate", ylim=c(0,1),
          # xlab = experimentLabels,
          main = main,
          las=2 , cex.axis = 1,
          col = testTable.colors
  )
  ## Draw horizontal grid
  abline(h=seq(from = 0, to = 1, by = 0.1), col="darkgrey", lty="dotted")
  abline(h=seq(from = 0, to = 1, by= 0.05 ), lty=2)
  meanPermlabels <- apply(testing.error.rates[permTestExperiments], 1, mean)

  ## Mean misclassification rate for all the label-permuted tests
  abline(h= mean(meanPermlabels), col="red", lwd=3 , lty= 3)

  ## Draw a line with the theoretically computed random expected misclassification rate.
  ## In principle, this should fit the mean of the permutation tests.
  ## TO BE ONE BUT WE NEED TO EITHER PASS THE DATASET OR ASSOCIATE THIS METHOD TO THE CLASS.
  # abline(h = , col="red", lwd=3 , lty= 3)

  ## Plot legend
  legend("topleft", lwd = 4,
         bty="o", bg = "white",
         legend = c("Actual data", "Permuted Labels"),
         col = c(expColor, permColor),
         cex = 0.8, pch =0.2)
  par(mar=save.margins)

  if (!is.null(boxplotFile)){
    silence <- dev.off()
  }
}

