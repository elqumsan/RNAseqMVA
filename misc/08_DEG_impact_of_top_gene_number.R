
####################################################################################################
####################################################################################################
#### Experiment: measure hit rates with increasing number of variables ordered by DEG p-value ####
###################################################################################################

##### define the file to store memory Image for " the Number of DEG Ordered" test #####
image.dir <- file.path( project.parameters$global$dir$memoryImages, parameters$recountID)
dir.create(image.dir, showWarnings = FALSE, recursive = TRUE)
image.file <- file.path(image.dir, paste(sep = "", "train_test_no._of_DEG_ordered_",parameters$recountID , ".Rdata"))


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
        if(dataSetName == "DEGdataSets"){
          var.numbers <- c( 10, 40, 80, 200, 400, 1000,
                            seq(from=3000, to=nrow(studyCases[[recountID]]$datasetsForTest$DEGdataSets$DESeq2$orderedDataTable )-1, by = 70000), nrow(studyCases[[recountID]]$datasetsForTest$DEGdataSets$DESeq2$orderedDataTable))

          studyCases[[recountID]]$parameters$var.numbers <- var.numbers
          studyCases[[recountID]]$datasetsForTest[[dataSetName]]$parameters <- parameters
        } # end if dataset is included in DEGdataset
      }
      #  print (studyCases[[recountID]]$parameters$dir$tablesDetail)
    }
  }
}

## Default for quick test without iterating over all cases
permute <- FALSE

#DEG.object <- DataTableWithDEG(studyCases$filtered)

if (project.parameters$global$compute) {
  message.with.time("Train/test all computations with constant training proportion :")

  train.test.results.all.DEG.ordered.per.classifier <- list()

  ##### Loop over recountID ######
  #   for (recountID in selectedRecountIDs) {
  parameters <- studyCases[[recountID]]$parameters
  DEGdataset <- studyCases[[recountID]]$datasetsForTest$DEGdataSets


  for (classifier in project.parameters$global$classifiers) {

  train.test.results.DEG <- list()

 #### Associate all computation with permuted and not permuted class lables ####
  for (permute in project.parameters$global$permute) {



    for (dataSetName in c("DESeq2","edgeR")) {

      DEG <- DEGdataset[[dataSetName]]
    #  DEG <- DEG.object$DEGdataSetsType[[dataset]]

    #  v  <- 1
      for(v in 1:length(studyCases[[recountID]]$parameters$var.numbers)){
        varnb <- studyCases[[recountID]]$parameters$var.numbers[v]

        ## For the time being we do this experiment only with log2 normalised counts
        ## since we saw that it improves the result with all variables
        # data.type <- paste(dataset$dataType, dataset$method, sep = "_")
        # data.table <- na.omit( as.data.frame(get(data.type)[["orderedDataTable"]]))
       if (dataSetName == "DESeq2"){
         selected.DEG.names <- rownames( DEGdataset$orderedDataTable[1:varnb,])

       } else{
         selected.DEG.names <- rownames( DEG$orderedDataTable[1:varnb,])

       }

        ## Make sure that we select gene names present in the selected data type
        ## (some genes may be filtered out or technical reasons)
        valid.DEG.names <- selected.DEG.names[selected.DEG.names %in% rownames(DEG$orderedDataTable)]

        DEG.object <- DEG$orderedDataTable[valid.DEG.names,]
         DEGdataset$dataTable <- DEG.object

        #DEG$DEG.object <- DEG.object
        #DataTableWithTrainTestSets(DEG$DEG.object )


        ## dim(counts)

        #### Define experiment prefix ####

        variable.type <- paste(sep="_", "DEG-",dataset, "top", varnb, "var")
        exp.prefix <- paste(sep="_", project.parameters$global$classifiers, studyCases[[recountID]]$ID , variable.type)
        if (permute) {
          exp.prefix <- paste(sep="_", exp.prefix,project.parameters$global$perm.prefix)
        }
        message (format(Sys.time(), "%Y-%m-%d_%H%M%S"), "\t", paste("Experiment prefix: ", exp.prefix, "DEG-Method", dataset))

        train.test.results.DEG[[exp.prefix]] <-
          IterateTrainingTesting(DEGdataset,
                         classifier=classifier,
                      #   iterations = parameters$iterations,
                          # data.type=data.type,
                         # classes=classes,
                        # trainIndex =get( paste(sep = ".", "DEG", deg.method))[["trainIndex"]], testIndex = get( paste(sep = ".","DEG", deg.method))[["testIndex"]],
                         #variable.type = variable.type,
                         #trainingProportion = parameters$trainingProportion,
                         permute = permute
                      #   file.prefix = exp.prefix
                      )
      } # end of for loop over nb.variables
    } # end  of for loop over DEG.methods
  } # end of for loop over permuted lables

  #### Print the results of the effect of the number of DEG ordered on the efficiancy of each classifier ####
  ErrorRateBoxPlot(experimentList = train.test.results.DEG,
                   classifier=classifier,
                 #  data.type = studyCases[[recountID]]$datasetsForTest$DEGdataSets$DESeq2$dataType,
                   data.type =  DEG$dataType,
                 #  variable.type = project.parameters$global$variables.type[2],
                   main = paste(  "impact of number of variables sorted according DEG", "\n",
                                studyCases[[recountID]]$ID,";",
                                project.parameters$global$iteration,
                                  "iterations;",
                                  "DEG.data.type;" ,
                                  "ordered_variables",sep = ""))

  train.test.results.all.DEG.ordered.per.classifier[[recountID]][[classifier]] <- train.test.results.DEG

    } # end loop over classifiers

  #### Save an image of the results to enable reloading them withouht recomputing everything ####
  if (project.parameters$global$save.image) {
    save.image(file = image.file)
  }

  ##### if compution not required, you can load the image file without any computations ####

#  } # end iteration over recountID
} # end of compute condition

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

## For each experiment, the results were stored in an element of the list train.test.results.
## Test the list of experiments done.
names(train.test.results.DEG)

#################################################################################################################
#### results for the imapct of the number of variables (genes) sorted according to DEG ####
#################################################################################################################
# ErrorRateBoxPlot(experimentList = train.test.results.DEG,
#                  classifier=classifier,
#
#                  variable.type = "ordered_variables_2",
#                  main = paste(  "impact of number of variables sorted according DEG", "\n",
#                                 parameters$recountID,";",
#                                 parameters$iterations,
#                                 "iterations;",
#                                 "DEG.data.type;" ,
#                                 "ordered_variables",sep = ""))

