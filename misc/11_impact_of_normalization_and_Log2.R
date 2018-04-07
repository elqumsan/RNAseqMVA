#' @title study the effect of Log2 onto the efficiency of the classifier
#' @author Mustafa ABUELQUMSAN and Jacques van Helden
#' @description this script measures the impact of Log2 to improve the effeciancy
#' of classification. We start by calculateing the miscalssification errors when the classifier
#' is applied to the raw counts directly without any prior treatment. We then perform some
#' "Log2-transformed" treatment onto these raw count and study the miscalssification errors
#' in order to evaluate if normalization improves the accuracy.
#'
#' @param countTable
#' @param classes
#' @param classfier
#' @param trainingproportion
#' @param k # number of neighbour for the classfier
#'
#'

rawTable <- loaded$originalExperiment$countTable
stat.raw <- list()

# M = data.frame(matrix(rnorm(100000),nrow=500))
##### some basic statistics to exhibit the nature of the raw data #####
stat.raw$min = apply(rawTable,2,min)
stat.raw$fQuan = apply(rawTable,2,quantile,prob=0.25)
stat.raw$median = apply(rawTable,2,median)
stat.raw$tQaun = apply(rawTable,2,quantile,prob=0.75)
stat.raw$max = apply(rawTable,2,max)

stat.raw$mean = apply(rawTable,2,mean)
stat.raw$sd = apply(rawTable,2,sd)

# par(mfrow=c(1,2))

message.with.time("\t\tDrawing plots describing original count table statistics")

file.prefix <- paste("main_stat.raw_",parameters$recountID,sep = "")
boxplot.file <- file.path(dir.NormImpact, paste(file.prefix,"boxplot.pdf", sep = "_"))
pdf(file = boxplot.file)
boxplot(cbind(stat.raw$min, stat.raw$fQuan, stat.raw$median, stat.raw$tQaun, stat.raw$max, stat.raw$mean, stat.raw$sd) ,
        ylim= c(0,50000), names=c("min","F.Qaunt","median","T.Qaunt","max","mean","sd"), las=1, cex.axis = 0.7,
        xlab= c(paste("Sammary statistics for the raw Table",parameters$recountID)))
silence <- dev.off()

file.prefix <- paste("diff_bet_stat.raw_", parameters$recountID, sep = "")
boxplot.file <- file.path(dir.NormImpact, paste(file.prefix, "boxplot.pdf", sep = "_"))
pdf(file = boxplot.file)
boxplot(cbind(stat.raw$fQuan - stat.raw$min, stat.raw$median - stat.raw$min, stat.raw$tQaun - stat.raw$min, stat.raw$max - stat.raw$tQaun) , ylim= c(0,50000),
        names= c("fQaun-min", "median-min", "tQuan-min", "max-tQuan"), cex.axis = 0.7,
        xlab= c(paste("difference between some statistics\n","in Raw table" ,parameters$recountID)))
silence <- dev.off()

##### cheking the No. libsum for each classes in the Raw Table  ##################

# data.frame(libsum=apply(loaded$countTable, 1, sum), class=loaded$classes)

# x <- data.frame(q3 = apply(rawCounts$Counts, 1, quantile, q=0.75), sum = apply(rawCounts$Counts, 1, sum), class=loaded$classes)
#
# x <- data.frame(q3 = apply(rawCounts$Counts, 1, quantile, q=0.75), sum = apply(rawCounts$Counts, 1, sum), class=loaded$classes)

libsum=apply(loaded$originalExperiment$countTable, 2, sum)
x <- data.frame(libsum=apply(loaded$originalExperiment$countTable, 2, sum), class=loaded$originalExperiment$classLabels)

head(x)


file.prefix <- paste("libsum.raw_",parameters$recountID,sep = "")
boxplot.file <- file.path(dir.NormImpact, paste(file.prefix,"boxplot.pdf", sep = "_"))
pdf(file = boxplot.file)
save.margins <- par("mar")
par(mar = c(14,5,5,1))
boxplot(libsum ~ class , data=x, horizontal=TRUE, las=1, cex.axis=0.7,
        xlab= c(paste("Libsum statistics for the raw Table\n","distribution for each class in",parameters$recountID)))
par(mar=save.margins)
silence <- dev.off()


##### some basic statistics to exhibit the nature of the log2norm data #####
log2normTable <- loaded$log2norm$counts
stat.log2norm <- list()

stat.log2norm$min = apply(log2normTable,2,min)
stat.log2norm$fQuan = apply(log2normTable,2,quantile,prob=0.25)
stat.log2norm$median = apply(log2normTable,2,median)
stat.log2norm$tQaun = apply(log2normTable,2,quantile,prob=0.75)
stat.log2norm$max = apply(log2normTable,2,max)

stat.log2norm$mean = apply(log2normTable,2,mean)
stat.log2norm$sd = apply(log2normTable,2,sd)

# par(mfrow=c(1,2))
message.with.time("Drawing plots describing log2norm table statistics")

file.prefix <- paste("main_stat.log2norm_",parameters$recountID,sep = "")
boxplot.file <- file.path(dir.NormImpact, paste(file.prefix,"boxplot.pdf", sep = "_"))
pdf(file = boxplot.file)
boxplot(cbind(stat.log2norm$min, stat.log2norm$fQuan, stat.log2norm$median, stat.log2norm$tQaun, stat.log2norm$max, stat.log2norm$mean, stat.log2norm$sd) ,
        ylim= c(0,100), names=c("min","F.Qaunt","median","T.Qaunt","max","mean","sd"), las=1, cex.axis = 0.7,
        xlab= c(paste("Sammary statistics for the log2norm Table",parameters$recountID)))
silence <- dev.off()

file.prefix <- paste("diff_bet_stat.log2norm_", parameters$recountID, sep = "")
boxplot.file <- file.path(dir.NormImpact, paste(file.prefix, "boxplot.pdf", sep = "_"))
pdf(file = boxplot.file)
boxplot(cbind(stat.log2norm$fQuan - stat.log2norm$min, stat.log2norm$median - stat.log2norm$min, stat.log2norm$tQaun - stat.log2norm$min, stat.log2norm$max - stat.log2norm$tQaun) , ylim= c(0,100),
        names= c("fQaun-min", "median-min", "tQuan-min", "max-tQuan"), cex.axis = 0.7,
        xlab= c(paste("difference between some statistics\n","in log2norm table" ,parameters$recountID)))
silence <- dev.off()

##### cheking the No. libsum for each classes in the log2norm Table  ##################

# data.frame(libsum=apply(loaded$countTable, 1, sum), class=loaded$classes)

# x <- data.frame(q3 = apply(log2norm$Counts, 1, quantile, q=0.75), sum = apply(log2norm$Counts, 1, sum), class=loaded$classes)
#
# x <- data.frame(q3 = apply(rawCounts1, 1, quantile, q=0.75), sum = apply(rawCounts1, 1, sum), class=loaded$classes)

x <- data.frame(libsum=apply(loaded$log2norm$counts, 2, sum), class=loaded$log2norm$classLabels)

head(x)

#

file.prefix <- paste("libsum.log2norm_",parameters$recountID,sep = "")
boxplot.file <- file.path(dir.NormImpact, paste(file.prefix,"boxplot.pdf", sep = "_"))
pdf(file = boxplot.file)
save.margins <- par("mar")
par(mar = c(14,5,5,1))
boxplot(libsum ~ class , data=x, horizontal=TRUE, las=1, cex.axis=0.7,
        xlab= c(paste("Libsum statistics for the log2norm Table\n","distribution for each class in",parameters$recountID)))
par(mar=save.margins)
silence <- dev.off()


#############################################################################################
###### Draw some plots to showing effects of some treatments(Norm and log2) on the Count Table ####
#############################################################################################
message.with.time("\tDrawing plots describing count table statistics and distribution")

#### Compute a trimmed mean: suppress the 5% top and bottom values ####
if (parameters$compute) {
  message.with.time("Computing trimmed mean of normalized counts")
  x <- unlist(loaded$norm$counts)
  q0.05 <- quantile(x = x, probs = 0.05, na.rm=TRUE)
  q0.95 <- quantile(x = x, probs = 0.95, na.rm=TRUE)
  trimmed <- (x[x > q0.05 & x < q0.95])
  suppressed.proportion <- 1 - length(trimmed)/length(x)
} else {
  message.with.time("Skipping a trimmed mean")
}





## histogram of normalized counts. Zoom on the representative part of the histogram
file.prefix <- file.path(dir.NormImpact, paste(parameters$recountID, "_counts_norm_perc75_hist.pdf", sep = ""))

pdf(file=file.prefix,
  width = 8, height = 8)
hist(unlist(loaded$norm$counts), breaks=10, las=1,
     xlab="normalized counts", ylab="Occurrences", col="grey",
     main=paste(parameters$recountID, "Normalised count distrib"))
abline(v=mean(unlist(loaded$norm$counts)), lwd=1, col="darkgreen") # mark the mean
abline(v=median(unlist(loaded$norm$counts)), lwd=1, col="blue") # mark the median
abline(v=mean(unlist(trimmed)), lwd=1, col="purple") # mark the trimmed mean
legend("topright",lwd=2,
       legend=c("mean", "median", "trimmed mean"),
       col=c("darkgreen", "blue", "purple"))
silence <- dev.off()


#### Compute a trimmed mean: suppress the 5% top and bottom values ####
if (parameters$compute) {
  message.with.time("Computing trimmed mean of log2normalized counts")
  x <- unlist(loaded$log2norm$counts)
  q0.05 <- quantile(x = x, probs = 0.05, na.rm=TRUE)
  q0.95 <- quantile(x = x, probs = 0.95, na.rm=TRUE)
  log2.trimmed <- (x[x > q0.05 & x < q0.95])
  suppressed.proportion <- 1 - length(log2.trimmed)/length(x)
} else {
  message.with.time("Skipping a trimmed mean")
}

file.prefix <- file.path(dir.NormImpact, paste(parameters$recountID, "_counts_norm_perc75_and_log2_hist.pdf", sep = ""))

pdf(file=file.prefix,
    width = 8, height = 8)
hist(unlist(loaded$log2norm$counts), breaks=100, las=1,
     xlab="log2 Normalized counts", ylab="Occurrences", col="grey",
     main=paste(parameters$recountID, " log 2 Normalised count distribution"))
abline(v=mean(unlist(loaded$log2norm$counts)), lwd=1, col="darkgreen") # mark the mean
abline(v=median(unlist(loaded$log2norm$counts)), lwd=1, col="blue") # mark the median
abline(v=mean(unlist(log2.trimmed)), lwd=1, col="purple") # mark the trimmed mean
legend("topright",lwd=2,
       legend=c("mean", "median", "trimmed mean"),
       col=c("darkgreen", "blue", "purple"))
silence <- dev.off()


######  interpretation the difference between the mean and the third quartile raw count #####
file.prefix <- file.path(dir.NormImpact, paste(parameters$recountID, "_rawcounts_mean_vs_Q3.pdf", sep = ""))
pdf(file= file.prefix,
  width = 8, height = 8)
plot(x=apply(loaded$originalExperiment$countTable, 2, mean),
     y=signif(digits=3, apply(loaded$originalExperiment$countTable, 2, quantile, 0.75)),
     main="raw counts: Percentile 75 versus mean",
     xlab="Mean counts per sample",
     ylab="Percentile 75",
     col=sampleColors,
     las=1,
     panel.first=grid())
legend("topleft", legend=names(classColors), col=classColors, pch=1, cex=0.8)
silence <- dev.off()

######  interpretation the difference between the mean and the third quartile log2counts #####
file.prefix <- file.path(dir.NormImpact, paste(parameters$recountID, "_log2norm_mean_vs_Q3.pdf", sep = ""))
pdf(file= file.prefix,
    width = 8, height = 8)
plot(x=apply(loaded$log2norm$counts, 1, mean),
     y=signif(digits=3, apply(loaded$log2norm$counts, 1, quantile, 0.75)),
     main="log2norm counts: Percentile 75 versus mean",
     xlab="Mean counts per sample",
     ylab="Percentile 75",
     col=sampleColors,
     las=1,
     panel.first=grid())
legend("topleft", legend=names(classColors), col=classColors, pch=1, cex=0.8)
silence <- dev.off()

message.with.time(" finishing from drawing plots describing count table statistics and distribution")

# ## Running parameters
# parameters <- list(
#   iterations = 10, ## Number of iterations for the classiifers
#   reload = TRUE
# )

# if (parameters$reload == TRUE) {
#   ################################################################################
#   ## Save an image of the memory, so I can reload it later to avoid re-running all the analyses
#   message("Loading memory image ")
#   load(file = image.file)
# }


# loading required libraries
# requiredCRAN <- c('class', "randomForest","broom", "devtools")
# RequiredCRANPackages(requiredCRAN)
#
#
# requiredBioconductor <- c("recount")
# RequiredBioconductorPackages(requiredBioconductor)
################################################################
## Load a count Table from recount-experiment, merge counts per sample
## and apply some pre-filtering (remove zero-variance and near-zero-variance genes).
# recountID <- "SRP048759"

# message("Loading count table from recount", "; recountID = ", recountID)
# loaded <- loadCounts(recountID = recountID, mergeRuns = T, classColumn = "tissue")
# rawCounts <- loaded$countTable ## Note: one row per sample, one column per gene
# dim(rawCounts)

################################################################
## Assign a specific color to each sammple according to its class
# pheno <- loaded$phenoTable
# classes <- loaded$classes
# distinct.classes <- as.vector(unique(loaded$classes))


##########################################################################################
############ the analysis KNN classifier with the Log2-normalised Real Dataset
###########################################################################################

# message("KNN classifier with log2-normalised real data set , ", parameters$iterations, " iterations.")
#
# iteratedKnnTestingErrorRatesLog2Norm <- data.frame()
# for (i in 1:parameters$iterations) {
#   # computing the testing errors rate for the KNN classfier
#   oneTestKnn <- MisclassificationEstimate(log2normCounts , loaded$classes ,trainingProportion = 2/3, classifier = "knn" )
#   iteratedKnnTestingErrorRatesLog2Norm <- rbind(iteratedKnnTestingErrorRatesLog2Norm, oneTestKnn$stats)
# }
#
# boxplot(iteratedKnnTestingErrorRatesLog2Norm[,"testing.error.rate"], horizontal = T, ylim= c(0,1), las=1)
#

###########################################################################################################
##########  Testing KNN classifier performance with the DESeq2 and edgeR ordaring variables for
##########  Log2-transformed-Normalized Real dataset
###########################################################################################################
## we thus will using the ordered variables which are the most significance resulting from DEG,
## concering with padj.

#message("KNN classifier with DESeq2 and edgeR ordaring real data set , ", parameters$iterations, " iterations.")

#DEG.DESeq <-  DEGordering(loaded$countTable , loaded$classes, method = "DESeq2")
#DEG.edgeR  <- DEGordering(loaded$countTable, loaded$classes, method = "edgeR")

# sorted.log2.transformed.edgeR <-log2normCounts[, DEG.DESeq$geneOrder]
# sorted.log2.transformed.DESeq <- log2normCounts[, DEG.edgeR$geneOrder]
#
# iteratedKnn.test.edgeR.Log2Norm <- data.frame()
#
# nb.variables <- c(3, 5, 7, 10, 15, 20, 30, 50, 100, 200, 500, 1000, 2000, 5000,10000, 20000 )
#
# #iteratedVarnb <- vector()
# orderediteratedKnn.test.DESeq.log2Norm.by.varnb <- list()
# v  <- 1
# for(v in 1:length(nb.variables)){
#
#   varnb <- nb.variables[v]
#   var.label <- paste(sep="_", "DESeq2_top", varnb)
#   message("DESeq2-ordered ", v, "/", length(nb.variables), ": ",varnb," variables")
#
#   iteratedKnn.test.DESeq.Log2Norm <- data.frame()
#   for(i in 1:parameters$iterations){
#     #test.edgeR <-MisclassificationEstimate(countTable = sorted.log2.transformed.edgeR[, 1:varnb], classes = classes ,classifier = "knn" )
#     test.DESeq <-MisclassificationEstimate(countTable = sorted.log2.transformed.DESeq[, 1:varnb], classes = classes, classifier = "knn")
#     #iteratedKnn.test.edgeR.Log2Norm <- rbind(iteratedKnn.test.edgeR.Log2Norm, test.edgeR$stats)
#     iteratedKnn.test.DESeq.Log2Norm <- rbind(iteratedKnn.test.DESeq.Log2Norm, test.DESeq$stats)
#   } # end for iteration
#
#   # Store the full result table for each var number
#   orderediteratedKnn.test.DESeq.log2Norm.by.varnb[[var.label]] <- iteratedKnn.test.DESeq.Log2Norm
#
#   # Collect the error rates in a single table with one row per iteration and one col per number of variables
#   if (v == 1) {
#
#     orderediteratedKnn.test.DESeq.log2Norm.err <- data.frame(iteratedKnn.test.DESeq.Log2Norm$testing.error.rate)
#     colnames(orderediteratedKnn.test.DESeq.log2Norm.err) <- var.label
#   } else {
#     orderediteratedKnn.test.DESeq.log2Norm.err[, var.label] <-
#       iteratedKnn.test.DESeq.Log2Norm$testing.error.rate
#
#   }
#
#   # orderediteratedKnn.test.DESeq.log2Norm <- cbind( as.list(orderediteratedKnn.test.DESeq.log2Norm),
#   #                                     as.list( iteratedKnn.test.DESeq.Log2Norm[,"testing.error.rate"]),  )
#
#   #iteratedmeanOfEachIteration <- append( iteratedmeanOfEachIteration ,meanOfEachIteration )
#   #iteratedVarnb <- append(iteratedVarnb, varnb )
# } # end for number of variables
#
# orderediteratedKnn.test.DESeq.log2Norm.err[,paste(sep="", "All variables (", ncol(log2normCounts),")")] <-
#   iteratedKnnTestingErrorRatesLog2Norm$testing.error.rate
#
#
# boxplot(orderediteratedKnn.test.DESeq.log2Norm.err, horizontal=TRUE, las=1,
#         ylim=c(0,1))
# abline(v=seq(from=0, to=1, by=0.1))
#
# ################################################################################
# #### Ploting the Miscalssification error rate with log2-transformed
# write.table(file = file.path(table.dirs[classifier], paste(parameters$recountID, "Impact_of_Log2_transformed_onto_ordering_variables.tsv")),
#             x = orderediteratedKnn.test.DESeq.log2Norm.err )
#
# pdf(file = file.path(figure.dirs[classifier], paste(parameters$recountID, "Impact_of_Log2_transformed_onto_ordering_variables", ".pdf")),
#     width = 9, height = 7 )
# save.margins <- par("mar")
# par(mar = c(3, 10, 7, 1))
# iteratedKnnTestingErrorRatesAllLog2NormColors <- c(1:length(orderediteratedKnn.test.DESeq.log2Norm.err))
# names(iteratedKnnTestingErrorRatesAllLog2NormColors ) <- names(orderediteratedKnn.test.DESeq.log2Norm.err)
# iteratedKnnTestingErrorRatesAllLog2NormColors[grep(names(iteratedKnnTestingErrorRatesAllLog2NormColors),pattern = "top")] <-"#9900FF"
# iteratedKnnTestingErrorRatesAllLog2NormColors[grep(names(iteratedKnnTestingErrorRatesAllLog2NormColors),pattern = "ALL")] <- "grey"
# boxplot(orderediteratedKnn.test.DESeq.log2Norm.err, las=1, horizontal = T, ylim= c(0,1),
#         col= iteratedKnnTestingErrorRatesAllLog2NormColors, cex.axis= 0.8,
#         xlab="Misclassification Rate",
#         main =paste("Misclassification Rates\n", parameters$iterations, "iterations"))
# abline(v=seq(0 ,1 , 0.05), col= "#DDDDDD", lty= "dashed")
# abline(v=seq(0, 1, 0.01), col="#DDDDDD", lty="dashed")
# par(mar = save.margins)
#
# silence <- dev.off()
#
# ##########################################################################################
# ############ the analysis KNN classifier with the Log2-normalised Permuted Dataset
# ###########################################################################################
#
# message("KNN classifier with log2-normalised permuted data set , ", parameters$iterations, " iterations.")
#
# iteratedKnnTestingErrorRatesLog2Norm <- data.frame()
# for (i in 1:parameters$iterations) {
#   # computing the testing errors rate for the KNN classfier
#   oneTestKnn <- MisclassificationEstimate(log2normCounts , sample(loaded$classes) ,trainingProportion = 2/3, classifier = "knn" )
#   iteratedKnnTestingErrorRatesLog2NormPermuted <- rbind(iteratedKnnTestingErrorRatesLog2NormPermuted, oneTestKnn$stats)
# }
#
#
#
#
# ###########################################################################################################
# ##########  Testing KNN classifier performance with the DESeq2 and edgeR ordaring variables for
# ##########  Log2-transformed-Normalized Permuted dataset
# ###########################################################################################################
# ## we thus will using the permuted ordered variables which are the most significance resulting from DEG,
# ## concering with padj.
#
# message("KNN classifier with DESeq2 and edgeR ordaring Permuted data set , ", parameters$iterations, " iterations.")
#
# sorted.log2.transformed.edgeR <- log2normCounts[, DEG.edgeR$geneOrder]
# sorted.log2.transformed.DESeq <- log2normCounts[, DEG.DESeq$geneOrder]
#
# iteratedKnn.test.DESeq.Log2NormPermuted <- data.frame()
# iteratedKnn.test.DESeq.Log2NormPermuted.by.varnb <- list()
# nb.variables <- c(3, 5, 7, 10, 15, 20, 30, 50, 100, 200, 500, 1000, 2000, 5000,10000, 20000)
#
# v <- 1
# for(v in 1:length(nb.variables)){
#   varnb <- nb.variables[v]
#   var.label <- paste("DESeq2_top", varnb, sep = "_")
#   message("DESeq2-ordered- Permuted ", v, "/", length(nb.variables), ": ",varnb," variables")
#   for(i in parameters$iterations){
#     oneTestlog2Normalised <-MisclassificationEstimate(countTable = sorted.log2.transformed.DESeq[, 1:varnb], classes = sample(loaded$classes), classifier = "knn")
#     iteratedKnn.test.DESeq.Log2NormPermuted <- rbind(iteratedKnn.test.edgeR.Log2NormPermuted , oneTestlog2Normalised$stats)
#   }
#   iteratedKnn.test.DESeq.Log2NormPermuted.by.varnb[[var.label]] <- iteratedKnn.test.DESeq.Log2NormPermuted
#
#   if(v == 1){
#
#     orderediteratedKnn.testDESeq.log2Norm.error <- data.frame( iteratedKnn.test.DESeq.Log2NormPermuted$testing.error.rate)
#     colnames(orderediteratedKnn.testDESeq.log2Norm.error) <- var.label
#   }
#   else {
#     orderediteratedKnn.testDESeq.log2Norm.error[, var.label] <- iteratedKnn.test.DESeq.Log2NormPermuted$testing.error.rate
#       }
#
#
# }
#
#
# orderediteratedKnn.testDESeq.log2Norm.error[, paste(sep = "","All Variables (" ,ncol(log2normCounts), ")" )] <-
#   data.frame(iteratedKnnTestingErrorRatesLog2NormPermuted$testing.error.rate )
# ################################################################################
# #### Ploting the Miscalssification error rate with log2-transformed
# write.table(file = file.path(table.dirs[classifier], paste(parameters$recountID,"Impact_of_Log2_transformed_onto_ordering_variables_permuted.tsv")),
#             x = orderediteratedKnn.testDESeq.log2Norm.error )
# pdf(file = file.path(figure.dirs[classifier], paste("Impact_of_Log2_transformed_onto_ordering_variables_permuted", ".pdf")),
#     width = 9, height = 7)
# save.margins <- par("mar")
# par(mar = c(3, 10, 7, 1))
# iteratedKnnTestingErrorRatesLog2NormColors <- c(1:length(orderediteratedKnn.testDESeq.log2Norm.error))
# names(iteratedKnnTestingErrorRatesLog2NormColors) <- names(orderediteratedKnn.testDESeq.log2Norm.error)
# iteratedKnnTestingErrorRatesLog2NormColors[grep(names(iteratedKnnTestingErrorRatesLog2NormColors),pattern = "top")] <-"#9900FF"
# iteratedKnnTestingErrorRatesLog2NormColors[grep(names(iteratedKnnTestingErrorRatesLog2NormColors),pattern = "ALL")] <- "grey"
# boxplot(orderediteratedKnn.testDESeq.log2Norm.error, las=1, horizontal = T, ylim= c(0,1),
#         col= iteratedKnnTestingErrorRatesLog2NormColors, cex.axis= 0.8,
#         xlab="Misclassification Rate",
#         main =paste("Misclassification Rates\n", parameters$iterations, "iterations"))
# abline(v=seq(0 ,1 , 0.05), col= "#DDDDDD", lty= "dashed")
# abline(v=seq(0, 1, 0.01), col="#DDDDDD", lty="dashed")
# par(mar = save.margins)
# silence <- dev.off()
#
