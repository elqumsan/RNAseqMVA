####################################################################################################
####################################################################################################
#### experiment: measure hit rates with increasing number of variables ordered by imporatnce of variables by random forest ####
###################################################################################################


##### define the file to store memory Image for " the Number of PCs" test #####
image.dir <- file.path( parameters$dir$memoryImages, parameters$recountID)
dir.create(image.dir, showWarnings = FALSE, recursive = TRUE)
image.file <- file.path(image.dir, paste(sep = "", "train_test_no._of_VIs_",parameters$recountID , ".Rdata"))


## For debug: reset the parameteres for all the study cases
## This is used to re-run the analyses on each study case after having changed some parameters in the yaml-specific configuration file
reload.parameters <- TRUE
if (reload.parameters) {
  project.parameters <- yaml.load_file(configFile)
  project.parameters <- initParallelComputing(project.parameters)
  if(exists("studyCases")) {
    for (recountID in names(studyCases)) {
      parameters <- initRecountID(recountID, project.parameters)
      studyCases[[recountID]]$parameters <- parameters
      for (dataSetName in names(studyCases[[recountID]]$datasetsForTest)) {
        ViRf.numbers <- c( 10, 40, 80, 200, 400, 1000,
                        seq(from=3000, to=nrow(studyCases[[recountID]]$datasetsForTest$log2normViRf$sigviRf )-1, by = 70000), nrow(studyCases[[recountID]]$datasetsForTest$log2normViRf$sigviRf))

        studyCases[[recountID]]$parameters$ViRf.numbers <- ViRf.numbers
        studyCases[[recountID]]$datasetsForTest[[dataSetName]]$parameters <- parameters
      }
      #  print (studyCases[[recountID]]$parameters$dir$tablesDetail)
    }
  }
}

# ##### define the file to store memory Image for " the Number of DEG Ordered" test #####
# image.dir <- file.path( parameters$dir$memoryImages, parameters$recountID)
# dir.create(image.dir, showWarnings = FALSE, recursive = TRUE)
# image.file <- file.path(image.dir, paste(sep = "", "train_test_no._of_v.importance_ordered_",parameters$recountID , ".Rdata"))
#
#
# ##### Define all the train indices for all the iterations, in order to using the same training\testing parts with different classifiers and data type. #####
# if (parameters$identicalTrainTest) {
#   ## New option: define all the train indices for all the iterations, in order to use the same training/testing sets between dfferent classifiers and data types
#   trainIndices <- list()
#   for(i in 1:parameters$iterations) {
#     n <- nrow(rawCounts$Counts)
#     train.size <- round(parameters$trainingProportion * n)
#     trainIndices [[i]] <- sample(1:n, size = train.size, replace = FALSE)
#   }
# } else {
#   ## First option: select different indices at each experiment
#   trainIndices = NULL
# }
#
#
#
# # dim(counts)
# # View(counts)

## Default for quick test without iterating over all cases


if (project.parameters$global$compute) {

  message.with.time("Train/test all computations with constant training proportion :")
  train.test.results.v.importance.ordered.per.classifier <- list()
 ##### Loop over recountID ######
#   for (recountID in selectedRecountIDs) {

    parameters <- studyCases[[recountID]]$parameters
    dataset <- studyCases[[recountID]]$datasetsForTest$log2normViRf

  for (classifier in project.parameters$global$classifiers) {

    ## List to store all results
    train.test.results.importance.varaibles <- list()
    #    train.test.results.all.PCs <- list()

    message.with.time("Train/test all computations with constant training proportion :",
                      signif(project.parameters$global$trainingProportion, digits = 3) )

    message.with.time("\tTrain/test, ", "; classifier=", classifier)




    #### Associate each analysis of real data with a permutation test ####
    permute <- FALSE
    for (permute in project.parameters$global$permute) {

    #  v.importance <- get("ordered.dataTable.by.importance")
      data.type <- "log2normViRf"
      v  <- 1
      for(v in 1:length(studyCases[[recountID]]$parameters$ViRf.numbers)){
        varnb <- studyCases[[recountID]]$parameters$ViRf.numbers[v]

        ## For the time being we do this experiment only with log2 normalised counts
        ## since we saw that it improves the result with all variables
       # data.type <- parameters$data.types["V.importance"]
       # data.table <- na.omit( as.data.frame(get(data.type)[["orderedDataTable"]]))
        selected.v.importance <- dataset$sigviRf[1:varnb,]
        selected.v.importance.names <-rownames(selected.v.importance)
        ## Make sure that we select gene names present in the selected data type
        ## (some genes may be filtered out or technical reasons)
        valid.v.importance.names <- selected.v.importance.names[selected.v.importance.names %in% rownames(dataset)]

        dataset.from.ViRf <- dataset$sigviRf[valid.v.importance.names,]

        #### Run classifier with the most importance variables (raw counts) ####

      ###### Define experiment prefix #######
        variable.type <- paste(sep = "_","v.import.rf","top", varnb,"var")
        exp.prefix <-
          paste(sep = "_", project.parameters$global$classifiers, studyCases[[recountID]]$ID , variable.type)
        if (permute) {
          exp.prefix <- paste(sep = "_", exp.prefix,project.parameters$global$perm.prefix)
        }# end if permuted class

        message (format(Sys.time(), "%Y-%m-%d_%H%M%S"), "\t", "Experiment prefix: ", exp.prefix)

        train.test.results.importance.varaibles[[exp.prefix]] <-
          IterateTrainingTesting (
            dataset,
           # classes = classes,
           # trainIndices = trainIndices,
            # trainIndex = sample(log2norm$trainIndex),
            # testIndex = sample(log2norm$testIndex),
            #data.type = parameters$data.types["V.importance"],
            classifier = classifier,
            #variable.type = variable.type,
          #  trainingProportion = parameters$trainingProportion,
           # file.prefix = exp.prefix,
            permute = permute
            # k = parameters$knn$k,
            # verbose = parameters$verbose
          )


      } # end of for loop nb.varaibles importance
      #  } # end of iterative of Variable importance Numbers
    } # end for loop permutation



    #### Ploting the Misclassification Error rate for ordered varaibles by the importances of each variables #####  ####
    ErrorRateBoxPlot(experimentList = train.test.results.importance.varaibles,
                     classifier = classifier,
                     data.type = dataset$dataType,
                     main = paste(sep="",
                                  classifier, ": impact of number of variables sorted according\n\t ","the most importance variables by random forest,", "\n",
                                  parameters$recountID, ", ",
                                  parameters$iterations, " iterations, ","\n",
                                  dataset$dataType))

      train.test.results.v.importance.ordered.per.classifier[[recountID]][[classifier]] <- train.test.results.importance.varaibles

  } # end of loop over classifiers

  #### Save an image of the results to enable reloading them withouht recomputing everything ####
  if (parameters$save.image) {
    save.image(file = image.file)
  }

  ##### if compution not required, you can load the image file without any computations ####

 # } ## end over recountID
} # end if Computed statent


# else {
#   # message.with.time("Skipping train/test analysis with DEG-ordered top variables")
#   # reload previous results if exist
#   if (file.exists(image.file)) {
#     message ("Reloading memory image ", image.file)
#     load(image.file)
#   } else {
#     stop("Cannot reload memory image file ", image.file)
#   }
#
# }

###############################################################################################
#### What is better to using ordered varaibles by the importances of each variables  ####
# ErrorRateBoxPlot(experimentList = train.test.results.all.importance.varaibles,
#                  classifier = classifier,
#                  data.type = parameters$data.types["V.importance"],
#                  main = paste(sep="",
#                               classifier, ": all importance variables,", "\n",
#                               parameters$recountID, ", ",
#                               parameters$iterations, " iterations, ","\n",
#                               data.type = "diverse data type"))

