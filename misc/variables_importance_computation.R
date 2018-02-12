########### Computation the importace of all vaiables in raw count table and ordered it by the importance ########
countTable <-  rawCounts$Counts
classes <- loaded$classes

# n <- nrow(rawCounts$Counts) ## Number of observations (samples)
# trainingProportion <- 0.66
# train.size <- round(n * parameters$trainingProportion)
#
# ## Random selection of indices for the training set
# trainIndex <- sort(sample(1:n, size=train.size))
# ## Use remaining indices for the testing set
# testIndex <- setdiff(1:n, trainIndex)

rf.model  <- randomForest(
  x = countTable,
  y =  as.factor( classes),
  xtest = countTable, importance = T, keep.forest = T)

variable.importance <- importance(rf.model,type = 1,scale = F)
#variable.importance[,1]

ordered.varaible.importance <-order(variable.importance[,1],decreasing = T)

ordered.countTable.by.importace <-countTable[, ordered.varaible.importance]

sig.variables <- round(ncol(ordered.countTable.by.importace) * 0.75)
ordered.countTable.by.importance  <- ordered.countTable.by.importace[, 1:sig.variables]


