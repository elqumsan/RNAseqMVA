
################################################################
## QUESTION: What is imapct of number of PCAs-transformed data?,
# which is the better to use a subsets of the first components or all the components ?

## Choice of the classifier
classifier <- "svm"

data.type <- "log2norm.prcomp.centred"

## iterate over permutation status
pc.numbers <- c(2, 3, 4, 5, 6, 7,
                seq(from=10, to=ncol(counts)-1, by = 10), ncol(counts))
pc.nb <- 4 ## Default or quick test


if (parameters$compute) {
  ## List to store all results
  train.test.results.No.PCs <- list()

  message.with.time("Train/test all computations with constant training proportion :",
                    signif(parameters$trainingProportion, digits = 3) )
  #for (k in parameters$knn$k) {
  #    for (classifier in classifier) {
  message.with.time("\tTrain/test, k=", parameters$knn$k, "; classifier=", classifier)
  ## iterate over data types
  #   for (data.type in parameters$data.types) {
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
    #    } else if (data.type == data.type){
    #      counts <- counts
  } else {
    stop(data.type, " is not a valid type, Supported: raw, norm, log2, and *pcr*")
  } # end else of other data.type

  ## iterate over permutation status
  pc.numbers <- c(2, 3, 4, 5, 6, 7,
                  seq(from=10, to=ncol(counts)-1, by = 10), ncol(counts))
  pc.nb <- 4 ## Default or quick test



  ## Associate each analysis of real data with a permutation test
  for (permute in c(FALSE, TRUE)) {

    ## Run classifier withb all variables (log2-transformed log counts)
    # exp.prefix <-
    #   paste(sep = "_", classifier, "log2norm", "allvars", "K", parameters$knn$k)
    # if (permute) {
    #   exp.prefix <- paste(sep = "_", exp.prefix, perm.prefix)
    # }# end if permuted class
    #

    # train.test.results[[exp.prefix]] <-
    #   one.experiment (
    #     countTable = log2normCounts,
    #     classes = classes,
    #     data.type = "log2normCounts",
    #     classifier = classifier,
    #     #variable.type = variable.type,
    #     trainingProportion = parameters$trainingProportion,
    #     file.prefix = exp.prefix,
    #     permute = permute,
    #     k = parameters$knn$k,
    #     verbose = parameters$verbose
    #   )


    ## Iterate over PC numbers
    for (pc.nb in pc.numbers) {

      ## Select the first N principal components
      first.pcs <- data.frame(counts[,1:pc.nb])
      colnames(first.pcs) <- colnames(counts)[1:pc.nb]
      variable.type <- paste(sep="", "PC_1-", pc.nb)

      ## define experiment prefix
      exp.prefix <-
        paste(sep = "_", classifier, data.type, "nb_of_PCs", pc.nb)
      if (permute) {
        exp.prefix <- paste(sep = "_", exp.prefix, perm.prefix)
      }# end if permuted class

      train.test.results.No.PCs[[exp.prefix]] <-
        one.experiment (
          countTable = first.pcs,
          classes = classes,
          data.type = data.type,
          classifier = classifier,
          variable.type = "PC_comp",
          trainingProportion = parameters$trainingProportion,
          file.prefix = exp.prefix,
          permute = permute,
          k = parameters$knn$k,
          verbose = parameters$verbose
        )

    } # end of iterative of PC Numbers
  } # end of permutation
} # end of computation

###############################################################################################
## Print the results of the effect of the number of PCs on the efficiancy of KNN classifier
ErrorRateBoxPlot(experimentList = train.test.results.No.PCs,
                 classifier = classifier,
                 variable.type = "number_of_PCs",
                 main = paste("Impact of the number of PCs on", classifier, "\n",
                              parameters$iterations , "iterations,",
                              data.type))

