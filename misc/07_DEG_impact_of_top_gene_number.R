
####################################################################################################
####################################################################################################
#### Second experiment: measure hit rates with increasing number of variables ordered by DEG p-value ####
###################################################################################################

## Choice of the classifier

classifier <- "knn"


## Choice of the coutns
#data.type <- "log2norm.prcomp.centred"
#data.type <- "log2norm"

# dim(counts)
# View(counts)

## Default for quick test without iterating over all cases
permute <- FALSE


if (parameters$compute) {
  message.with.time("Starting classification")

  train.test.results.DEG <- list()

 #### Associate all computation with permuted and not permuted class lables ####
  for (permute in c(FALSE, TRUE)) {

    for (deg.method in parameters$deg.methods) {
      DEG <- get(paste(sep="", "DEG.", deg.method))

      v  <- 5
      for(v in 1:length(parameters$nb.variables)){
        varnb <- parameters$nb.variables[v]

        ## For the time being we do this experiment only with log2 normalised counts
        ## since we saw that it improves the result with all variables
        data.type <- paste("DEG", deg.method, sep = ".")
        data.table <- na.omit( as.data.frame(get(data.type)[["orderedCountTable"]]))
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
                         trainIndex =get( paste(sep = ".", "DEG", deg.method))[["trainIndex"]], testIndex = get( paste(sep = ".","DEG", deg.method))[["testIndex"]],
                         variable.type = variable.type,
                         trainingProportion = parameters$trainingProportion,
                         permute = permute,
                         file.prefix = exp.prefix)
      }
    }
  }
  #  } # end for the assiciate the analysis for each classifier
}  else {
  message.with.time("Skipping train/test analysis with DEG-ordered top variables")
}

## For each experiment, the results were stored in an element of the list train.test.results.
## Test the list of experiments done.
names(train.test.results.DEG)

#################################################################################################################
#### results for the imapct of the number of variables (genes) sorted according to DEG ####
#################################################################################################################
ErrorRateBoxPlot(experimentList = train.test.results.DEG,
                 classifier=classifier,

                 variable.type = "ordered_variables_2",
                 main = paste(  "impact of number of variables sorted according DEG", "\n",
                                parameters$recountID,";",
                                parameters$iterations,
                                "iterations;",
                                "DEG.data.type;" ,
                                "ordered_variables",sep = ""))

