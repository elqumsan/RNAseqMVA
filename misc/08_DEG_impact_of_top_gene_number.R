
####################################################################################################
####################################################################################################
#### Second experiment: measure hit rates with increasing number of variables ordered by DEG p-value ####
###################################################################################################

##### define the file to store memory Image for " the Number of DEG Ordered" test #####
image.dir <- file.path( parameters$dir$memoryImages, parameters$recountID)
dir.create(image.dir, showWarnings = FALSE, recursive = TRUE)
image.file <- file.path(image.dir, paste(sep = "", "train_test_no._of_DEG_ordered_",parameters$recountID , ".Rdata"))


# ##### Define all the train indices for all the iterations, in order to using the same training\testing parts with different classifiers and data type. #####
# if (parameters$identicalTrainTest) {
#   ## New option: define all the train indices for all the iterations, in order to use the same training/testing sets between dfferent classifiers and data types
#   trainIndices <- list()
#   for(i in 1:parameters$iterations) {
#     n <- nrow(log2norm$Counts)
#     train.size <- round(parameters$trainingProportion * n)
#     trainIndices [[i]] <- sample(1:n, size = train.size, replace = FALSE)
#   }
# } else {
#   ## First option: select different indices at each experiment
#   trainIndices = NULL
# }
## Choice of the coutns
#data.type <- "log2norm.prcomp.centred"
#data.type <- "log2norm"

# dim(counts)
# View(counts)

## Default for quick test without iterating over all cases
permute <- FALSE

DEG.object <- countTableWithDEG(loaded$filtered)

if (parameters$compute) {
  message.with.time("\t\tStarting classification")

  train.test.results.all.DEG.ordered.per.classifier <- list()

  for (classifier in parameters$classifiers) {

  train.test.results.DEG <- list()

 #### Associate all computation with permuted and not permuted class lables ####
  for (permute in c(FALSE, TRUE)) {


    for (dataset in DEG.object$DEG.datasets) {
      DEG <- dataset

      v  <- 5
      for(v in 1:length(parameters$nb.variables)){
        varnb <- parameters$nb.variables[v]

        ## For the time being we do this experiment only with log2 normalised counts
        ## since we saw that it improves the result with all variables
        # data.type <- paste(dataset$dataType, dataset$method, sep = "_")
        # data.table <- na.omit( as.data.frame(get(data.type)[["orderedCountTable"]]))
        dataset$
        selected.DEG.names <- DEG$geneOrder[1:varnb]

        ## Make sure that we select gene names present in the selected data type
        ## (some genes may be filtered out or technical reasons)
        valid.DEG.names <- selected.DEG.names[selected.DEG.names %in% colnames(data.table)]

        counts <- data.table[,valid.DEG.names]


        ## dim(counts)

        #### Define experiment prefix ####
        variable.type <- paste(sep="_", "DEG", deg.method, "top", varnb)
        exp.prefix <- paste(sep="_", classifier,  parameters$recountID, variable.type)
        if (permute) {
          exp.prefix <- paste(sep="_", exp.prefix, perm.prefix)
        }
        message (format(Sys.time(), "%Y-%m-%d_%H%M%S"), "\t", "Experiment prefix: ", exp.prefix)

        train.test.results.DEG[[exp.prefix]] <-
          one.experiment(countTable=counts,
                         data.type=data.type,
                         classifier=classifier,
                         classes=classes,
                        # trainIndex =get( paste(sep = ".", "DEG", deg.method))[["trainIndex"]], testIndex = get( paste(sep = ".","DEG", deg.method))[["testIndex"]],
                         variable.type = variable.type,
                         trainingProportion = parameters$trainingProportion,
                         permute = permute,
                         file.prefix = exp.prefix)
      } # end of for loop over nb.variables
    } # end  of for loop over DEG.methods
  } # end of for loop over permuted lables

  #### Print the results of the effect of the number of DEG ordered on the efficiancy of each classifier ####
  ErrorRateBoxPlot(experimentList = train.test.results.DEG,
                   classifier=classifier,
                   data.type = parameters$data.types["DEG"],
                   variable.type = "ordered_variables_2",
                   main = paste(  "impact of number of variables sorted according DEG", "\n",
                                  parameters$recountID,";",
                                  parameters$iterations,
                                  "iterations;",
                                  "DEG.data.type;" ,
                                  "ordered_variables",sep = ""))

  train.test.results.all.DEG.ordered.per.classifier[[classifier]] <- train.test.results.DEG

    } # end loop over classifiers

  #### Save an image of the results to enable reloading them withouht recomputing everything ####
  if (parameters$save.image) {
    save.image(file = image.file)
  }

  ##### if compution not required, you can load the image file without any computations ####


}  else {
  # message.with.time("Skipping train/test analysis with DEG-ordered top variables")
  # reload previous results if exist
  if (file.exists(image.file)) {
    message ("Reloading memory image ", image.file)
    load(image.file)
  } else {
    stop("Cannot reload memory image file ", image.file)
  }

}

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

