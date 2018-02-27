#################### Load counts and pheno ####################
## Load a count Table from recount-experiment, merge counts per sample
## and apply some pre-filtering (remove zero-variance and near-zero-variance genes).
if (parameters$compute) {
  rawCounts <- list()
  message.with.time("Loading count table from recount", "; recountID = ", parameters$recountID)
  loaded <- loadCounts(recountID = parameters$recountID, mergeRuns = TRUE ,
                       classColumn = parameters$classColumn,
                       minSamplesPerClass = parameters$minSamplesPerClass)
  rawCounts1 <- loaded$countTable ## Note: one row per sample, one column per gene
  # dim(rawCounts)

  ## Assign a specific color to each sammple according to its class
  pheno <- loaded$phenoTable
  classes <- loaded$classes
  geo.characteristics <- loaded$geo.characteristics
  # table(classes)
  # length(classes)
  distinct.classes <- as.vector(unique(loaded$classes))
  write.table(rawCounts, file = paste(tsv.dir,"/rawCounts_",parameters$recountID,".tsv", sep = ""), row.names = FALSE, sep = "\t")
  characteristics.string <- unlist(lapply(pheno$characteristics, paste, collapse="; "))
  pheno.data.frame <- data.frame(pheno$sample, characteristics.string)

  dim(pheno.data.frame)
  # View(pheno.data.frame)
  write.table(pheno.data.frame, file = paste(tsv.dir, "/pheno_",parameters$recountID,".tsv",sep = ""),
              row.names = FALSE, sep = "\t")
} else {
  message.with.time("Skipping data loading")
}

# Check the dimensions of the count table

dim(rawCounts1)
rawCounts$Counts <- na.omit(rawCounts1)
if(nrow(rawCounts$Counts) != length(classes)){
  stop(" the Number of samples should be eqaul to Number of classes")
}

######### sptiting the rawCounts dataset for the train set and test set #########
n <- nrow(rawCounts$Counts) ## Number of observations (samples)
train.size <- round(n * parameters$trainingProportion)

## Random selection of indices for the training set
trainIndex <- sort(sample(1:n, size=train.size))
## Use remaining indices for the testing set
testIndex <- setdiff(1:n, trainIndex)

rawCounts$trainIndex <- trainIndex
rawCounts$testIndex  <- testIndex

## Number of samples per class
print(loaded$samples.per.class)

## I start by assigning one systematic color(number) to each class,
## To make sure that each class has a color even with different datasets analsed in the efuture.
classColors <- 1:length(distinct.classes)
names(classColors) <- distinct.classes

##JvH: Mustafa, these colors are specific for one dataset.
## In a future phase we will discuss about how to manage dataset-specific parameters,
## but not before your seminar Let us just keep this in mind for later.

## Then Assign some specific colors which evoke particular sampletypes
if (parameters$recountID == "SRP048759") {
  classColors["Heparinised.blood"] <- "#BB0000"
  classColors["Bone.marrow"] <- "#4488FF"
} else if (parameters$recountID == "SRP042620") {
  ## TO DO : define colors for the multi-group breast cancer dataset
  classColors["Breast.Cancer.Cell.Line"] <- "red"
  classColors["ER..Breast.Cancer.Primary.Tumor"] <- "darkblue"
  classColors["Triple.Negative.Breast.Cancer.Primary.Tumor"] <- "purple"
  classColors["Uninvolved.Breast.Tissue.Adjacent.to.ER..Primary.Tumor"] <- "green"
  classColors["Uninvolved.Breast.Tissue.Adjacent.to.TNBC.Primary.Tumor"] <- "black"
}
# print(classColors)
## Assign colors per sample according to their class
sampleColors <- classColors[loaded$classes]
names(sampleColors) <- rownames(loaded$phenoTable)
# print(sampleColors)

##### Normalize the counts without log2 transformation (first test) #####
##
## Note: this method takes a table with one column per sample and one
## row per gene, we thus have to transpose the raw count table.
if (parameters$compute) {
  message.with.time("Normalizing counts based on 75th percentile")
  norm <- NormalizeCounts(t(loaded$countTable), method = "quantile", quantile=0.75, log2 = FALSE)
  normCounts <- t(norm$normCounts) ## Transpose normalized counts for classifiers
  dim(normCounts)

} else {
  message.with.time("Skipping normalisation")
}

#### Export normalized counts ####
dir.create(dir.NormImpact, recursive = TRUE, showWarnings = FALSE)
list.files(dir.NormImpact)
normCounts.file <- file.path(dir.NormImpact, paste(sep="", parameters$recountID, "_normalized_counts.tsv"))
message.with.time("Exporting normalised counts to file ", "\n", normCounts.file)
if (parameters$save.tables) {
  write.table(sep="\t", quote=FALSE, row.names = TRUE, col.names=NA,
              x = round(t(normCounts), digits=2),
              file = normCounts.file)
  write.table(x = round(t(normCounts), digits = 3), file = paste(tsv.dir,"/NormCounts_",parameters$recountID,".tsv", sep = ""),
              row.names = FALSE, sep = "\t")
} else {
  message.with.time("Skipping saving of normalized counts table")
}

##### Normalize counts with log2 transformation (second test) #####
##
## Note: this method takes a table with one column per sample and one
## row per gene, we thus have to transpose the raw count table.
log2norm <- list()
if (parameters$compute) {
  message.with.time("Normalizing counts based on 75th percentile + log2 transformation")
  log2normCounts <- NormalizeCounts(t(loaded$countTable), method = "quantile", quantile=0.75,
                              log2 = TRUE, epsilon=0.1)
  Counts <- na.omit( as.data.frame(t(log2normCounts$normCounts)))
  log2norm$Counts <- Counts
  dim(log2norm$Counts)

  if(nrow(log2norm$Counts) != length(classes)){
    stop(" the Number of samples in log2norm counts should be the same length of classes")
  }

  ######### sptiting the log2norm dataset for the train set and test set #########
  n <- nrow(log2norm$Counts) ## Number of observations (samples)
  train.size <- round(n * parameters$trainingProportion)

  ## Random selection of indices for the training set
  trainIndex <- sort(sample(1:n, size=train.size))
  ## Use remaining indices for the testing set
  testIndex <- setdiff(1:n, trainIndex)

  log2norm$trainIndex <- trainIndex
  log2norm$testIndex  <- testIndex

} else {
  message.with.time("Skipping normalisation with log2 transformation")
}

#### Export log2-transformed normalized counts ####
dir.create(dir.NormImpact, recursive = TRUE, showWarnings = FALSE)
log2normCounts.file <- file.path(dir.NormImpact, paste(sep="", parameters$recountID, "_log2_normalized_counts.tsv"))
message.with.time("Exporting log2 normalised counts to file ", "\n", log2normCounts.file)
if (parameters$save.tables) {
  write.table(sep="\t", quote=FALSE, row.names = TRUE, col.names=NA,
              x = round(digits=3, t(log2norm$Counts)),
              file = log2normCounts.file)
}

#### Compute a trimmed mean: suppress the 5% top and bottom values ####
if (parameters$compute) {
  message.with.time("Computing trimmed mean of normalized counts")
  x <- unlist(normCounts)
  q0.05 <- quantile(x = x, probs = 0.05, na.rm=TRUE)
  q0.95 <- quantile(x = x, probs = 0.95, na.rm=TRUE)
  trimmed <- (x[x > q0.05 & x < q0.95])
  suppressed.proportion <- 1 - length(trimmed)/length(x)
} else {
  message.with.time("Skipping a trimmed mean")
}


#### Compute a trimmed mean: suppress the 5% top and bottom values ####
if (parameters$compute) {
  message.with.time("Computing trimmed mean of log2normalized counts")
  x <- unlist(log2norm$Counts)
  q0.05 <- quantile(x = x, probs = 0.05, na.rm=TRUE)
  q0.95 <- quantile(x = x, probs = 0.95, na.rm=TRUE)
  log2.trimmed <- (x[x > q0.05 & x < q0.95])
  suppressed.proportion <- 1 - length(log2.trimmed)/length(x)
} else {
  message.with.time("Skipping a trimmed mean")
}


## Indicate that this script has finished running
message.with.time("finished executing 02_load_and_normalise_counts.R")

##### plotting some figures to explore the nuture of recount data set #####
message.with.time(" plotting some figures to explore distribution for the recount data set ",parameters$recountID)
source("misc/11_impact_of_normalization_and_Log2.R")

##### Exhibiting the geo charactiristics for the current project #####
message.with.time("Exhibit the geo charactiristics for such experiment: ", parameters$recountID, " in order to know the class lable
                  for such experiment")
head( geo.characteristics)
geo.characteristics.file <- file.path("~/RNAseqMVA_workspace", "data", parameters$recountID, "geo.characteristics.tsv")
write.table( geo.characteristics, file = geo.characteristics.file, quote = FALSE,
             row.names = FALSE, sep = "\t" )
