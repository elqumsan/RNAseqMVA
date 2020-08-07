#' @title Tune the parameters of the different classifiers for a given study case
#' @author Jacques van Helden
#' @param studyCase a studyCase object
#' @param dataType="TMM_log2" data type to use for the tuning. Must be one of the data types included in the 'datasetsForTest' attribute of the studyCase object
#' @param tuneSVM=TRUE if TRUE, tune parameters for SVM
#' @param svmCost=1 cost parameter values for SVM
#' @param svmDegree=1:4 degree values for SVM polynomial kernel
#' @param svmGamma=4^(-2:4) gamma values for SVM kernels (ignored for linear kernel)
#' @param tuneRandomForest=TRUE if TRUE, tune parameters for Random Forest
#' @param rfNtree=c(100,200,300,500) values to test for RF numbers of trees
#' @param rfMtry=c(50,100,200,300,500) values to test for RF mtry parameter
#' @param tuneKNN=TRUE if TRUE, tune parameters for KNN
#' @param knnK=c(1,2,3,4,5,7,10) number of neighbours for KNN
#' @param knnL=0 minimum vote for definite decision for KNN
#' @param plotResults=FALSE if TRUE, run plot()  on the tuned objects
#' @import e1071
#' @import randomForest
#' @return a list of parameters cloned from the studyCase, added with the optimal parameters identified here
#' @export
TuneClassifiers <- function(studyCase,
                            dataType="TMM_log2_PC",
                            tuneSVM=TRUE,
#                            svmCost=4^(-2:2),
                            svmCost=1,
                            svmDegree=1:4,
                            svmGamma=4^(-2:4),
                            svmCoef0=4^(-2:4),
                            tuneRF=TRUE,
                            rfNtree=c(100,200,300,500),
                            rfMtry=c(50,100,200,300,500),
                            tuneKNN=TRUE,
                            knnK=c(1,2,3,4,5,7,10),
                            knnL=0,
                            plotResults=FALSE) {

  ## Prepare the result object
  result <- list()
  result$parameters <- studyCase$parameters
  result$tuneResults <- list()


  #### Data selection ####
  dataSet <- studyCase$datasetsForTest[[dataType]]
  if (is.null(dataSet)) {
    stop("dataType not supported for this studyCase: ", dataType,
         "\n\nSupported data types: ", paste(collapse = ", ", names(studyCase$datasetsForTest)))
  }
  names(dataSet)



  #### Tuning of the SVM parameters ####
  if (tuneSVM) {
    message("\tTuneClassifiers()\tTuning SVM")
    result$tuneResults$svm <- list()

    #### Parameters for the linear kernel ####
    message.with.time("\t\tTuning SVM with linear kernel\t", studyCase$ID, " ", studyCase$parameters$feature)
    linearTuning <- tune.svm(x = t(dataSet$dataTable),
                             y = as.factor(dataSet$classLabels),
                             cost = svmCost,
                             kernel = "linear")
    result$tuneResults$svm$linearTuning <- linearTuning
    if (plotResults) {
      plot(linearTuning, main = "SVM tuning - linear kernel")
    }


    #### Parameters for the polynomial kernel ####
    message.with.time("\t\tTuning SVM with polynomial kernel\t", studyCase$ID, " ", studyCase$parameters$feature)
    polynomialTuning <- tune.svm(x = t(dataSet$dataTable),
                                 y = as.factor(dataSet$classLabels),
                                 cost = svmCost,
                                 gamma = svmGamma,
                                 degree = svmDegree,
                                 coef0 = svmCoef0,
                                 kernel = "polynomial")
    result$tuneResults$svm$polynomialTuning <- polynomialTuning
    if (plotResults) {
      plot(polynomialTuning, main = "SVM tuning - polynomial kernel")
    }


    #### Parameters for the radial kernel ####
    message.with.time("\t\tTuning SVM with radial kernel\t", studyCase$ID, " ", studyCase$parameters$feature)
    radialTuning <- tune.svm(x = t(dataSet$dataTable),
                             y = as.factor(dataSet$classLabels),
                             cost = svmCost,
                             gamma = svmGamma,
                             kernel = "radial")
    result$tuneResults$svm$radialTuning <- radialTuning
    if (plotResults) {
      plot(radialTuning, main = "SVM tuning - radial kernel")
    }

    #### Parameters for the sigmoid kernel ####
    message.with.time("\t\tTuning SVM with sigmoid kernel\t", studyCase$ID, " ", studyCase$parameters$feature)
    sigmoidTuning <- tune.svm(x = t(dataSet$dataTable),
                              y = as.factor(dataSet$classLabels),
                              cost = svmCost,
                              gamma = svmGamma,
                              coef0 = svmCoef0,
                              kernel = "sigmoid")
    result$tuneResults$svm$sigmoidTuning <- sigmoidTuning
    if (plotResults) {
      plot(sigmoidTuning, main = "SVM tuning - sigmoid kernel")
    }

  }

  #### Tune Random Forest parameters ####
  if (tuneRandomForest) {
    message("\tTuneClassifiers()\tTuning RF")
    result$tuneResults$rf <- list()

    ## Tune random forest with e1071::tune.randomForest() function
    tunedRandomForest  <-
      tune.randomForest(x = t(dataSet$dataTable),
                        y = as.factor(dataSet$classLabels),
                        nodesize = 1,
                        mtry = rfMtry,
                        ntree = rfNtree)
    result$tuneResults$rf$tunedRandomForest <- tunedRandomForest
    if (plotResults) {
      plot(tunedRandomForest, main = "RF tuning")
    }


    # ## use randomForest::tuneRF to tune RF
    # tunedRF <- randomForest::tuneRF(
    #   x = t(dataSet$dataTable),
    #   y = as.factor(dataSet$classLabels),
    #   #                     mtryStart = 8,
    #   ntreeTry = 500,
    #   stepFactor = 1.5,
    #   improve = 0.025,
    #   trace = TRUE,
    #   plot = TRUE,
    #   doBest = FALSE)
    # result$tuneResults$rf$tunedRF <- tunedRF
    # if (plotResults) {
    #   plot(tunedRF, main = "RF tuning")
    # }

  }

  #### Tune KNN parameters ####
  if (tuneKNN) {
    message("\tTuneClassifiers()\tTuning KNN")
    result$tuneResults$knn <- list()

    tunedKNN <- tune.knn(x = t(dataSet$dataTable),
                         y = as.factor(dataSet$classLabels),
                         k = knnK,
                         l = knnL)

    result$tuneResults$knn$tunedKNN <- tunedKNN
    if (plotResults) {
#      plot(tunedKNN, main = "KNN tuning")
      plot(tunedKNN$performances$k,
           tunedKNN$performances$error,
           type = "b",
           main = paste(sep = "",
                        studyCase$ID,
                        " " , studyCase$parameters$feature,
                        " - ", studyCase$parameters$short_label,
                        "\n", "KNN tuning"),
           xlab = "k",
           ylab = "MER",
           las = 1,
           panel.first = grid(),
           col = "red"
           )
    }

  }

  return(result)

}

