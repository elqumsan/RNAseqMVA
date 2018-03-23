################################################################
#### QUESTION: What is imapct of number of PCAs-transformed data?, ####
# which is the better to use a subsets of the first components or all the components ?

##### define the file to store memory Image for " the Number of PCs" test #####
image.dir <- file.path( parameters$dir$memoryImages, parameters$recountID)
dir.create(image.dir, showWarnings = FALSE, recursive = TRUE)
image.file <- file.path(image.dir, paste(sep = "", "train_test_no._of_PCs_",parameters$recountID , ".Rdata"))

#data.type <- "log2norm.prcomp.centred"
#data.type <- "log2norm"
#data.type <- parameters$data.types["prcomp"]
data.type <- "log2norm.prcomp.centred.scaled"

##### Define all the train indices for all the iterations, in order to using the same training\testing parts with different classifiers and data type. #####
if (parameters$identicalTrainTest) {
  ## New option: define all the train indices for all the iterations, in order to use the same training/testing sets between dfferent classifiers and data types
  trainIndices <- list()
  for(i in 1:parameters$iterations) {
    n <- nrow(log2norm$Counts)
    train.size <- round(parameters$trainingProportion * n)
    trainIndices [[i]] <- sample(1:n, size = train.size, replace = FALSE)
  }
} else {
  ## First option: select different indices at each experiment
  trainIndices = NULL
}

#### iterate over permutation status ####
pc.numbers <- c(2, 3, 4, 5, 6, 7,
                seq(from=10, to=ncol(log2norm.prcomp.centred.scaled$x)-1, by = 10), ncol(log2norm.prcomp.centred.scaled$x))
pc.nb <- 4 ## Default or quick test


if (parameters$compute) {

  train.test.results.all.PCs.per.classifier <- list()

  for (classifier in parameters$classifiers) {

  ## List to store all results
  train.test.results.No.PCs <- list()

  message.with.time("Train/test all computations with constant training proportion :",
                    signif(parameters$trainingProportion, digits = 3) )
  #for (k in parameters$knn$k) {
  #    for (classifier in classifier) {
  message.with.time("\tTrain/test, k=", parameters$knn$k, "; classifier=", classifier)
  ## iterate over data types
  #   for (data.type in parameters$data.types) {
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
  #   #    } else if (data.type == data.type){
  #   #      counts <- counts
  # } else {
  #   stop(data.type, " is not a valid type, Supported: raw, norm, log2, and *pcr*")
  # } # end else of other data.type

  #### iterate over permutation status ####
  # pc.numbers <- c(2, 3, 4, 5, 6, 7,
  #                 seq(from=10, to=ncol(counts)-1, by = 10), ncol(counts))
  # pc.nb <- 4 ## Default or quick test



  #### Associate each analysis of real data with a permutation test ####
  for (permute in c(FALSE, TRUE)) {

    datasets <- list(loaded$filtered, loaded$norm ,loaded$log2norm)

    for (dataset in datasets) {

      dataset <- PCsWithTrainTestSets(dataset)
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


    #### Iterate over PC numbers ####
    for (pc.nb in dataset$pc.numbers) {

      ## Select the first N principal components
      first.pcs <- data.frame(dataset$PCs[,1:pc.nb])
      colnames(first.pcs) <- colnames(dataset$PCs)[1:pc.nb]
      dataset$first.pcs <- first.pcs
     # variable.type <- paste(sep="", "PC_1-", pc.nb)

      ## define experiment prefix
      exp.prefix <-
        paste(sep = "_", classifier, dataset$ID  , dataset$dataType, "nb_of_PCs", pc.nb)
      if (permute) {
        exp.prefix <- paste(sep = "_", exp.prefix, perm.prefix)
      }# end if permuted class

      train.test.results.No.PCs[[exp.prefix]] <-
        one.experiment (
          dataset,
        # classes = classes,
        #  trainIndex = log2norm.prcomp.centred.scaled$trainIndex, testIndex = log2norm.prcomp.centred.scaled$testIndex,
        #  data.type = data.type,
          classifier = classifier,
         # variable.type = "PC_comp",
          #trainingProportion = parameters$trainingProportion,
          file.prefix = exp.prefix,
          permute = permute,
          k = parameters$knn$k,
          verbose = parameters$verbose
        )
      } # end iteration over nb of PCs
    } # end of iteration over datasets
   } # end of permutation

  #### Print the results of the effect of the number of PCs on the efficiancy of each classifier classifier ####
  ErrorRateBoxPlot(experimentList = train.test.results.No.PCs,
                   classifier = classifier,
                   data.type = parameters$data.types["prcomp"],
                   #variable.type = "number_of_PCs",
                   main = paste("Impact of the number of PCs on,", classifier, "\n,",parameters$recountID,";",
                                parameters$iterations , "iterations,",
                                dataset$dataType , sep = ""),
                   variablesType = dataset$variablesType)

  train.test.results.all.PCs.per.classifier[[classifier]] <- train.test.results.No.PCs

   }  # end of loop over classifiers
 }  # end of if computation
#   #### Save an image of the results to enable reloading them withouht recomputing everything ####
#   if (parameters$save.image) {
#     save.image(file = image.file)
#   }
#
#   ##### if compution not required, you can load the image file without any computations ####
#
# } else {
#   # reload previous results if exist
#   if (file.exists(image.file)) {
#     message ("Reloading memory image ", image.file)
#     load(image.file)
#   } else {
#     stop("Cannot reload memory image file ", image.file)
#   }
#
# } # end else if loading image file
#
# # end of computation

###############################################################################################
#### Print the results of the effect of the number of PCs on the efficiancy of KNN classifier ####
# ErrorRateBoxPlot(experimentList = train.test.results.No.PCs,
#                  classifier = classifier,
#                  variable.type = "number_of_PCs",
#                  main = paste("Impact of the number of PCs on,", classifier, "\n,",parameters$recountID,";",
#                               parameters$iterations , "iterations,",
#                               data.type,sep = ""))

