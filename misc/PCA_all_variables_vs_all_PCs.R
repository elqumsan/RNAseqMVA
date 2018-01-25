
################################################################
## All variables versus all PCs.
##
## QUESTION: is it better to use the PCAs-transformed data, and, if so, is it better to use a subset of the first components or all the components ?
## For the time being we test this with only one classifier (KNN, default k)  but we will come back to it with other classifiers later.

## Choice of the classifier

classifier <- "rf"

## Choice of the coutns
#data.type <- "log2norm.prcomp.centred"
data.type <- "log2norm.prcomp.centred.scaled"

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
  if (data.type == "raw") {
    counts <- rawCounts
  } else if (data.type == "norm") {
    counts <- normCounts
  } else if (data.type == "log2") {
    counts <- log2normCounts
  } else if (grepl("prc", data.type)) {
    ## If we are not in the cases above, we assume that the data type is a
    ## PCA results, and we need to get the components.
    ##
    ## Create a variable name "counts" with the content of the variable whose
    ## name is given in "data.type"
    counts <- get(data.type)[["x"]]
  } else {
    stop(data.type, " is not a valid type, Supported: raw, norm, log2, and *pcr*")
  } # end else of other data.type


  ## Associate each analysis of real data with a permutation test
  for (permute in c(FALSE, TRUE)) {

    ## Run classifier withb all variables (log2-transformed log counts)
    exp.prefix <-
      paste(sep = "_", classifier, "log2norm", "allvars")
    if (permute) {
      exp.prefix <- paste(sep = "_", exp.prefix, perm.prefix)
    }# end if permuted class

    train.test.results.all.variables[[exp.prefix]] <-
      one.experiment (
        countTable = log2normCounts,
        classes = classes,
        data.type = "log2normCounts",
        classifier = classifier,
        #variable.type = variable.type,
        trainingProportion = parameters$trainingProportion,
        file.prefix = exp.prefix,
        permute = permute,
        k = parameters$knn$k,
        verbose = parameters$verbose
      )

    ## take all the principal components, and cast them to a data.frame
    first.pcs <- data.frame(counts)

    ## define experiment prefix
    exp.prefix <-
      paste(sep = "_", classifier, data.type)
    if (permute) {
      exp.prefix <- paste(sep = "_", exp.prefix, perm.prefix)
    }# end if permuted class

    train.test.results.all.variables[[exp.prefix]] <-
      one.experiment (
        countTable = first.pcs,
        classes = classes,
        data.type = data.type,
        classifier = classifier,
        variable.type = "all_PCs",
        trainingProportion = parameters$trainingProportion,
        file.prefix = exp.prefix,
        permute = permute,
        k = parameters$knn$k,
        verbose = parameters$verbose
      )

    #  } # end of iterative of PC Numbers
  } # end of permutation
} # end of computation

###############################################################################################
## What is better to using all PCs versus all variables with KNN classifier?
ErrorRateBoxPlot(experimentList = train.test.results.all.variables,
                 classifier = classifier,
                 main = paste(classifier, ": all variables vs all PCs", "\n",
                              parameters$iterations, "iterations,",
                              data.type = data.type))

