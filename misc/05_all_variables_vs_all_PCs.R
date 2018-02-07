
################################################################
##### All variables versus all PCs. #####
##
## QUESTION: is it better to use the PCAs-transformed data, and, if so, is it better to use a subset of the first components or all the components ?
## For the time being we test this with only one classifier (KNN, default k)  but we will come back to it with other classifiers later.
## IMPORTANT NOTE : i would like to pay your intention for here we should take " data.type, so that we will not give the user to
## choose the data.type in return we will pass the data.type for each experiment.
## Choice of the classifier

classifier <- "svm"

## Choice of the Counts
# data.type <- "log2norm.prcomp.centred"
# data.type <- "log2norm"


# dim(counts)
# View(counts)

## Default for quick test without iterating over all cases
permute <- FALSE

if (parameters$compute) {
  ## List to store all results
  train.test.results.all.variables <- list()
  train.test.results.all.PCs <- list()

  message.with.time("Train/test all computations with constant training proportion :",
                    signif(parameters$trainingProportion, digits = 3) )
  message.with.time("\tTrain/test, k=", parameters$knn$k, "; classifier=", classifier)

  ## Select the counts depending on the data type
  # if (data.type == "raw") {
  #   counts <- rawCounts
  # } else if (data.type == "norm") {
  #   counts <- normCounts
  # } else if (data.type == "log2norm") {
  #   counts <- log2norm
  # } else if (grepl("prc", data.type)) {
  #   ## If we are not in the cases above, we assume that the data type is a
  #   ## PCA results, and we need to get the components.
  #   ##
  #   ## Create a variable name "counts" with the content of the variable whose
  #   ## name is given in "data.type"
  #   counts <- get(data.type)[["x"]]
  #
  # } else {
  #   stop(data.type, " is not a valid type, Supported: raw, norm, log2, and *pcr*")
  # } # end else of other data.type


  #### Associate each analysis of real data with a permutation test ####
  for (permute in c(FALSE, TRUE)) {

    #### Run classifier with all variables (log2-transformed log counts) ####
    exp.prefix <-
      paste(sep = "_", classifier, parameters$recountID , parameters$data.type["log2norm"] , "allvars")
    if (permute) {
      exp.prefix <- paste(sep = "_", exp.prefix, perm.prefix)
    }# end if permuted class

    train.test.results.all.variables[[exp.prefix]] <-
      one.experiment (
        countTable = as.data.frame(log2norm$Counts),
        classes = classes,
        trainIndex = sample(log2norm$trainIndex),
        testIndex = sample(log2norm$testIndex),
        data.type = parameters$data.types["log2norm"],
        classifier = classifier,
        #variable.type = variable.type,
        trainingProportion = parameters$trainingProportion,
        file.prefix = exp.prefix,
        permute = permute,
        k = parameters$knn$k,
        verbose = parameters$verbose
      )

    #### take all the principal components, and cast them to a data.frame ####
    #first.pcs <- data.frame(counts)
     first.pcs <- get("log2norm.prcomp.centred.scaled")
    ## define experiment prefix
    exp.prefix <-
      paste(sep = "_", classifier, parameters$recountID , parameters$data.types["prcomp"])
    if (permute) {
      exp.prefix <- paste(sep = "_", exp.prefix, perm.prefix)
    }# end if permuted class

    train.test.results.all.variables[[exp.prefix]] <-
      one.experiment (
        countTable = first.pcs$x,
        classes = classes,
        trainIndex = sample( log2norm.prcomp.centred.scaled$trainIndex),
        testIndex = sample(log2norm.prcomp.centred.scaled$testIndex),
        data.type = parameters$data.types["prcomp"],
        classifier = classifier,
        variable.type = "all_PCs",
        trainingProportion = parameters$trainingProportion,

        file.prefix = exp.prefix,
        permute = permute,
        k = parameters$knn$k,
        verbose = parameters$verbose
      )


    #### take all the raw data without any nolmalization, and cast them to a data.frame ####
    ## we looking here to notice the ipmact of normalization into classifiers

    rawCounts1 <- na.omit(as.data.frame(get("rawCounts")[["Counts"]]))
    ## define experiment prefix
    exp.prefix <-
      paste(sep = "_", classifier, parameters$recountID , parameters$data.types["raw"])
    if (permute) {
      exp.prefix <- paste(sep = "_", exp.prefix, perm.prefix)
    }# end if permuted class

    train.test.results.all.variables[[exp.prefix]] <-
      one.experiment (
        countTable = rawCounts1,
        classes = classes,
        trainIndex= sample( rawCounts$trainIndex),
        testIndex = sample(rawCounts$testIndex),
        data.type = parameters$data.types["raw"],
        classifier = classifier,
        variable.type = "raw",
        trainingProportion = parameters$trainingProportion,

        file.prefix = exp.prefix,
        permute = permute,
        k = parameters$knn$k,
        verbose = parameters$verbose
      )

    #### take all the raw data ordered by DEG edgeR tool, and cast them to a data.frame ####
    ## we looking here to notice the ipmact of ordered variables by DEG edgeR tool into classifiers

    DEG.object <- get("DEG.edgeR")
    DEG.Counts <- na.omit(as.data.frame(DEG.object$orderedCountTable))
    sig.variables <- round(ncol(DEG.Counts) * 0.75)
    DEG.Counts <- DEG.Counts[, 1:sig.variables]
    ## define experiment prefix
    exp.prefix <-
      paste(sep = "_", classifier, parameters$recountID , parameters$data.types["DEG"])
    if (permute) {
      exp.prefix <- paste(sep = "_", exp.prefix, perm.prefix)
    }# end if permuted class

    train.test.results.all.variables[[exp.prefix]] <-
      one.experiment (
        countTable = DEG.Counts,
        classes = classes,
        trainIndex = sample( DEG.edgeR$trainIndex ) ,
        testIndex = sample( DEG.edgeR$testIndex),
        data.type = parameters$data.types["DEG"],
        classifier = classifier,
        variable.type = "DEG",
        trainingProportion = parameters$trainingProportion,
        file.prefix = exp.prefix,
        permute = permute,
        k = parameters$knn$k,
        verbose = parameters$verbose
      )

    #  } # end of iterative of PC Numbers
  } # end for loop permutation
} # end if statment computation

###############################################################################################
#### What is better to using all PCs versus all variables with KNN classifier? ####
ErrorRateBoxPlot(experimentList = train.test.results.all.variables,
                 classifier = classifier,
                 main = paste(sep="",
                              classifier, ": all variables vs all PCs,", "\n",
                              parameters$recountID, ", ",
                              parameters$iterations, " iterations, ","\n",
                              data.type = "diverse data type"))

