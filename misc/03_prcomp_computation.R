
#### Compute principal components for normalized log2 counts ####
message.with.time("Pre-processing by Principal Component analysis (PCA)")
log2norm.prcomp.centred <- prcomp( na.omit(log2norm$Counts), center = TRUE, scale. = FALSE)
log2norm.prcomp.centred.scaled <- prcomp(na.omit(log2norm$Counts), center = TRUE, scale. = TRUE)

######### sptiting the log2norm.prcomp.centres.scaled[["x"]] dataset for the train set and test set #########

n <- nrow(log2norm.prcomp.centred.scaled$x) ## Number of observations (samples)
train.size <- round(n * parameters$trainingProportion)

## Random selection of indices for the training set
trainIndex <- sort(sample(1:n, size=train.size))
## Use remaining indices for the testing set
testIndex <- setdiff(1:n, trainIndex)

log2norm.prcomp.centred.scaled$trainIndex <- trainIndex
log2norm.prcomp.centred.scaled$testIndex <- testIndex
## Indicate that this principal components analysis for log2 count has finished running
message.with.time("finished running Principal Component analysis (PCA) for normalized log2 counts")


