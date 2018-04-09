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
#' @param  permLegend legend for the one experiment with permuted class lables.
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
ErrorRateBoxPlot <- function(experimentList,
                             classifier,
                             data.type, # = "log2norm.prcomp.centred",
#                             variablesType,
                             main = paste(sep="", parameters$recountID, "; ", classifier,
                                          "\n", data.type,
                                          #"; ", variablesType,
                                          "; ", parameters$iterations, "iterations"),
                             expColor = "#00BBFF",
                             permColor = "grey",
                             expLegend = "Train/test",
                             permLegend = "Permuted labels",
                             ....) {


  # ## Define file to store the boxplot
#  output.dir <- figure.dirs[classifier]
  testing.error.rates <- data.frame()
  experimentNames <- names(experimentList)


  ## Make sure the output directory exists for storing the plots

  if (is.null(experimentNames) || length(experimentNames) < 1) {
    message.with.time("Not a single experiment for all variable composition")
  } else {
    exp <-1
    for (exp in 1:length(experimentNames)) {
      experiment <- experimentNames[exp]
      exp.result <- experimentList[[experiment]]
      if (exp == 1) {
        testing.error.rates <- data.frame(exp.result$testing.error.rate)
      } else {
        testing.error.rates <-
          cbind(testing.error.rates,
                exp.result$testing.error.rate)
      }
    } # end iterate the all variables experiment

    colnames(testing.error.rates) <- experimentNames
    rownames(testing.error.rates) <- 1:nrow(testing.error.rates)

    ## Define file to store the boxplot
    file.prefix <- paste(sep="_",
                         parameters$recountID,
                         classifier,
                         data.type,
#                         variablesType,
                         "iteration",
                         parameters$iterations)
    if (permute) {
      file.prefix <- paste(sep = "_", file.prefix, parameters$perm.prefix)
    }# end if permuted class
    boxplot.file <- file.path(parameters$dir$figures[classifier],
                              paste(sep = "", file.prefix, "_boxplot.pdf"))
    message.with.time("Boxplot file:\t", boxplot.file)
    pdf(file = boxplot.file, width = 3 + 0.2*length(experimentNames), height = 12)

    ## Define parameters for the boxplot
    save.margins <- par("mar")
    par(mar=c(14, 5, 5, 1))

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
           # xlab = experimentNames,
            main = main,
            las=2 , cex.axis = 0.7,
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
    legend("topright", lwd = 4,
           bty="o", bg = "white",
           legend = c("Actual data", "Permuted Labels"),
           col = c(expColor, permColor),
           cex = 1, pch =0.2)
    par(mar=save.margins)
    silence <- dev.off()
  } # end else all variables experiments
}

