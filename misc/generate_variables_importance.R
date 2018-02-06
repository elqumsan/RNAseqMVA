
rawCounts <- loaded$countTable
classes <- loaded$classes

n <- nrow(rawCounts) ## Number of observations (samples)
trainingProportion <- 0.66
train.size <- round(n * trainingProportion)

## Random selection of indices for the training set
trainIndex <- sort(sample(1:n, size=train.size))
## Use remaining indices for the testing set
testIndex <- setdiff(1:n, trainIndex)

rf.trained  <- randomForest(
  x = countTable[trainIndex, ],
  y =  as.factor( classes[trainIndex]),
  xtest = countTable[testIndex,], importance = T, keep.forest = T)

rf.trained




