#' @title Tune the parameters of the different classifiers for a given study case
#' @author Jacques van Helden
#' @param studyCase a studyCase object
#' @param dataType="TMM_log2" data type to use for the tuning. Must be one of the data types included in the 'datasetsForTest' attribute of the studyCase object
#' @param tuneSVM=TRUE if TRUE, tune parameters for SVM
#' @param tuneRF=TRUE if TRUE, tune parameters for Random Forest
#' @param tuneKNN=TRUE if TRUE, tune parameters for KNN
#' @param plotResults=TRUE if TRUE, run plot()  on the tuned objects
#' @import e1071
#' @import randomForest
#' @return a list of parameters cloned from the studyCase, added with the optimal parameters identified here
#' @export
TuneClassifiers <- function(studyCase,
                            dataType="TMM_log2",
                            tuneSVM=TRUE,
                            tuneRF=TRUE,
                            tuneKNN=TRUE) {

  ## Prepare the result object
  result <- list()
  result$parameters <- studyCase$parameters
  result$tuneResults <-


  #### Data selection ####
  dataSet <- studyCase$datasetsForTest[[dataType]]
  if (is.null(dataSet)) {
    stop("dataType not supported for this studyCase: ", dataType,
         "\n\nSupported data types: ", paste(collapse = ", ", names(studyCase$datasetsForTest)))
  }
  names(dataSet)



  #### Tuning of the SVM parameters ####
  if (tuneSVM) {

    #### Parameters for the linear kernel ####
    linearTuning <- tune.svm(x = t(dataSet$dataTable),
                             y = as.factor(dataSet$classLabels),
                             #                 cost = 4^(-2:2),
                             kernel = "linear")
    linearTuning
    plot(linearTuning, main = "SVM tuning - linear kernel")
    alarm()


    #### Parameters for the polynomial kernel ####
    polynomialTuning <- tune.svm(x = t(dataSet$dataTable),
                                 y = as.factor(dataSet$classLabels),
                                 #                 cost = 4^(-2:2),
                                 degree = 1:4,
                                 kernel = "polynomial")
    polynomialTuning
    plot(polynomialTuning, main = "SVM tuning - polynomial kernel")
    alarm()

    #### Parameters for the radial kernel ####
    radialTuning <- tune.svm(x = t(dataSet$dataTable),
                             y = as.factor(dataSet$classLabels),
                             #                 cost = 4^(-2:2),
                             gamma = 4^(-2:4),
                             coef0 = 4^(-2:4),
                             kernel = "radial")
    radialTuning
    plot(radialTuning, main = "SVM tuning - radial kernel")
    alarm()

    # obj3 <- best.svm(x = t(dataSet$dataTable),
    #                  y = as.factor(dataSet$classLabels), kernel = "linear")
    # obj3
    # plot(obj3)
  }

    #### Tune Random Forest parameters ####

  tunedRF <- randomForest::tuneRF(
    x = t(dataSet$dataTable),
    y = as.factor(dataSet$classLabels),
    #                     mtryStart = 8,
    ntreeTry = 500,
    stepFactor = 1.5,
    improve = 0.025,
    trace = TRUE,
    plot = TRUE,
    doBest = FALSE)
  plot(tunedRF)
  alarm()


}

