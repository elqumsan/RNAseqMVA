################################################################
#### QUESTION: What is imapct of number of PCAs-transformed data?, ####
# which is the better to use a subsets of the first components or all the components ?

##### define the file to store memory Image for " the Number of PCs" test #####
image.dir <- file.path( parameters$dir$memoryImages, parameters$recountID)
dir.create(image.dir, showWarnings = FALSE, recursive = TRUE)
image.file <- file.path(image.dir, paste(sep = "", "train_test_no._of_PCs_",parameters$recountID , ".Rdata"))


## For debug: reset the parameteres for all the study cases
## This is used to re-run the analyses on each study case after having changed some parameters in the yaml-specific configuration file
reload.parameters <- TRUE
if (reload.parameters) {
  project.parameters <- yaml.load_file(configFile)
  project.parameters <- initParallelComputing(project.parameters)
  if (exists("studyCases")) {
    for (recountID in names(studyCases)) {
      parameters <- initRecountID(recountID, project.parameters)
      studyCases[[recountID]]$parameters <- parameters
      for (dataSetName in names(studyCases[[recountID]]$datasetsForTest)) {
        pc.numbers <- c(2, 3, 4, 5, 6, 7,
                        seq(from=10, to=nrow(studyCases[[recountID]]$datasetsForTest$log2normPCs$dataTable)-1, by = 10), nrow(studyCases[[recountID]]$datasetsForTest$log2normPCs$dataTable))

        studyCases[[recountID]]$parameters$pc.numbers <- pc.numbers
        studyCases[[recountID]]$datasetsForTest[[dataSetName]]$parameters <- parameters
      }
      #  print (studyCases[[recountID]]$parameters$dir$tablesDetail)
    }
  }
}

if (project.parameters$global$compute) {

  train.test.results.all.PCs.per.classifier <- list()

  ## Loop over recountIDs
  for (recountID in selectedRecountIDs) {

    #  train.test.results.all.PCs.per.classifier[[recountID]] <- list() ## Instantiate an entry per recountID
    parameters <- studyCases[[recountID]]$parameters
    #studyCases$PCsVar <- PCsWithTrainTestSets(studyCases$filtered)
    dataset <- studyCases[[recountID]]$datasetsForTest$log2normPCs$dataTable

    for (classifier in project.parameters$global$classifiers) {

      ## List to store all results
      train.test.results.No.PCs <- list()

      message.with.time("Train/test all computations with constant training proportion :",
                        signif(project.parameters$global$trainingProportion, digits = 3) )
      #for (k in parameters$knn$k) {
      #    for (classifier in classifier) {
      message.with.time("\tTrain/test, ", "; classifier=", classifier)
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
      permute <- FALSE
      for (permute in project.parameters$global$permute) {
        data.type <- "log2normPCs"

        dataset   <-studyCases[[recountID]]$datasetsForTest[[data.type]]
        # datasets <- list(studyCases$filtered, studyCases$norm ,studyCases$log2norm)




        ## Run classifier withb all variables (log2-transformed log counts)
        # exp.prefix <-
        #   paste(sep = "_", classifier, "log2norm", "allvars", "K", parameters$knn$k)
        # if (permute) {
        #   exp.prefix <- paste(sep = "_", exp.prefix, perm.prefix)
        # }# end if permuted class
        #

        # train.test.results[[exp.prefix]] <-
        #   IterateTrainingTesting (
        #     dataTable = log2normCounts,
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
        for (pc.nb in studyCases[[recountID]]$parameters$pc.numbers) {

          ## Select the first N principal components
          first.pcs <- data.frame(t(dataset$prcomp$x[,1:pc.nb]))
          rownames(first.pcs) <- rownames(t(dataset$prcomp$x[,1:pc.nb]))
          dataset$dataTable <- first.pcs

          # dataset$firstPCs[[pc.nb]] <- first.pcs
          # variable.type <- paste(sep="", "PC_1-", pc.nb)

          ## define experiment prefix
          # currentOutputParameters <- outputParameters(dataset = dataset, classifier = classifier,permute = permute, createDir = TRUE)
          exp.prefix <-
            paste(sep = "_", classifier, dataset$ID  , dataset$dataType, "nb_of_PCs", pc.nb)
          if (permute) {
            exp.prefix <- paste(sep = "_", exp.prefix,parameters$perm.prefix)
          }# end if permuted class

          message (format(Sys.time(), "%Y-%m-%d_%H%M%S"), "\t", paste("Experiment prefix: ", exp.prefix))


          train.test.results.No.PCs[[exp.prefix]] <-
            IterateTrainingTesting (
              dataset,
              # classes = classes,
              #  trainIndex = log2norm.prcomp.centred.scaled$trainIndex, testIndex = log2norm.prcomp.centred.scaled$testIndex,
              #  data.type = data.type,
              classifier = classifier,
              # variable.type = "PC_comp",
              #trainingProportion = parameters$trainingProportion,
              # file.prefix = exp.prefix,
              permute = permute
              # k = parameters$knn$k,
              # verbose = parameters$verbose
            )
        } # end iteration over nb of PCs

      } # end of permutation


      #### Print the results of the effect of the number of PCs on the efficiancy of each classifier classifier ####
      ErrorRateBoxPlot(experimentList = train.test.results.No.PCs,
                       classifier = classifier,
                       #variable.type = "number_of_PCs",
                       main = paste("Impact of the number of PCs on,", classifier, "\n,",parameters$recountID,";",
                                    parameters$iterations , "iterations,",
                                    dataset$dataType , sep = "")
                       #variablesType = project.parameters$global$variables.type[[2]]
      )


      #  train.test.results.all.PCs.per.classifier[[recountID]][[classifier]] <- train.test.results.No.PCs
    }  # end of loop over classifiers

  } # end loop over recountIDs
}  # end of if computation

