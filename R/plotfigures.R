#' @title plot diagnostic figures for specific classes of objects
#' @author Mustafa AbuElQUmsan and Jacques van Helden
#' @description to explore the nature of raw data and showing the the distribution and main statistics for the raw data.
#' @param self is an object that is belonge to DataTableWithClasses.
#' @return utmost figures and visualisation about the main statistics and the distribution for the raw data.
#' @export
plotFigures  <- function(self,...){
  message("Exporting object of class", class(self), "to figures")
  UseMethod("plotFigures", self)


}


#' @title plot diagnostic figures for an object of the class DataTableWithClasses.
#' @author Mustafa AbuElQUmsan and Jacques van Helden
#' @description to explore the nature of raw data and showing the the distribution and main statistics for the raw data.
#' @param self is an object that is belonge to DataTableWithClasses.
#' @return utmost figures and visualisation about the main statistics and the distribution for the raw data.
#' @export
plotFigures.DataTableWithClasses <- function( self,
                                               plot.dir,
                                               file.prefix,
                                               dataType = self[["dataType"]],
                                               extension = ".pdf") {


  message( "Plotting the nature of the original, filtered and normalised count Data ", self[["ID"]], "for experiment" )
  message("\tfigure directory\t", plot.dir)
  message("\tFile prefix\t", file.prefix)


  rawTable <- self$dataTable
  # dim(rawTable)
  stat.raw <- list()
  # M = data.frame(matrix(rnorm(100000),nrow=500))
  ##### some basic statistics to exhibit the nature of the raw data #####
  stat.raw$min = apply(rawTable,1,min)
  stat.raw$fQuan = apply(rawTable,1,quantile,prob=0.25)
  stat.raw$median = apply(rawTable,1,median)
  stat.raw$tQaun = apply(rawTable,1,quantile,prob=0.75)
  stat.raw$max = apply(rawTable,1,max)

  stat.raw$mean = apply(rawTable,1,mean)
  stat.raw$sd = apply(rawTable,1,sd)

  # par(mfrow=c(1,2))

  message.with.time("\t\tDrawing plots describing original count table statistics")

  #file.prefix <- paste("main_stat.raw_",parameters$recountID,sep = "")
  boxplot.file <- file.path(plot.dir , paste(dataType ,"_",self[["ID"]], "_basic_stat", extension, sep = ""))
  pdf(file = boxplot.file)
  boxplot(cbind(stat.raw$min, stat.raw$fQuan, stat.raw$median, stat.raw$tQaun, stat.raw$max, stat.raw$mean, stat.raw$sd) ,
          ylim= c(0,100), names=c("min","F.Qaunt","median","T.Qaunt","max","mean","sd"), las=1, cex.axis = 0.7,
          xlab= c(paste("Sammary statistics for the raw Table",self[["ID"]])), ylab=" No. of genes")
  silence <- dev.off()

  file.prefix <- paste(dataType,"_", self[["ID"]], sep = "")
  boxplot.file <- file.path(plot.dir, paste(file.prefix, "differ_stat",extension, sep = "_"))
  pdf(file = boxplot.file)
  boxplot(cbind(stat.raw$fQuan - stat.raw$min, stat.raw$median - stat.raw$min, stat.raw$tQaun - stat.raw$min, stat.raw$max - stat.raw$tQaun) , ylim= c(0,50),
          names= c("fQaun-min", "median-min", "tQuan-min", "max-tQuan"), cex.axis = 0.7,
          xlab= c(paste("difference between some statistics\n","in Raw table" ,parameters$recountID)),
          ylab = "No. genes")
  silence <- dev.off()

  ##### cheking the No. libsum for each classes in the Raw Table  ##################

  # data.frame(libsum=apply(loaded$dataTable, 1, sum), class=loaded$classes)

  # x <- data.frame(q3 = apply(rawCounts$Counts, 1, quantile, q=0.75), sum = apply(rawCounts$Counts, 1, sum), class=loaded$classes)
  #
  # x <- data.frame(q3 = apply(rawCounts$Counts, 1, quantile, q=0.75), sum = apply(rawCounts$Counts, 1, sum), class=loaded$classes)

  libsum=apply(self[["dataTable"]], 1, sum)
  x <- data.frame(libsum=apply(self[["dataTable"]], 2, sum), class=self[["classLabels"]])

  head(x)


  file.prefix <- paste("libsum_",dataType,"_",self[["ID"]],sep = "")
  boxplot.file <- file.path(plot.dir, paste(file.prefix,extension, sep = ""))
  pdf(file = boxplot.file)
  save.margins <- par("mar")
  par(mar = c(14,5,5,1))
  boxplot(libsum ~ self[["classLabels"]] , data=x, horizontal=TRUE, las=1, cex.axis=0.7,
          xlab= c(paste("Libsum statistics for the raw Table\n","distribution for each class in",self[["ID"]])),
          ylab="Class Labels")
  par(mar=save.margins)
  silence <- dev.off()


  ##### some basic statistics to exhibit the nature of the log2norm data #####
  # log2normTable <- loaded$log2norm$counts
  # stat.log2norm <- list()
  #
  # stat.log2norm$min = apply(log2normTable,2,min)
  # stat.log2norm$fQuan = apply(log2normTable,2,quantile,prob=0.25)
  # stat.log2norm$median = apply(log2normTable,2,median)
  # stat.log2norm$tQaun = apply(log2normTable,2,quantile,prob=0.75)
  # stat.log2norm$max = apply(log2normTable,2,max)
  #
  # stat.log2norm$mean = apply(log2normTable,2,mean)
  # stat.log2norm$sd = apply(log2normTable,2,sd)
  #
  # # par(mfrow=c(1,2))
  # message.with.time("Drawing plots describing log2norm table statistics")
  #
  # file.prefix <- paste("main_stat.log2norm_",parameters$recountID,sep = "")
  # boxplot.file <- file.path(parameters$dir$NormalizationImpact, paste(file.prefix,"boxplot.pdf", sep = "_"))
  # pdf(file = boxplot.file)
  # boxplot(cbind(stat.log2norm$min, stat.log2norm$fQuan, stat.log2norm$median, stat.log2norm$tQaun, stat.log2norm$max, stat.log2norm$mean, stat.log2norm$sd) ,
  #         ylim= c(0,100), names=c("min","F.Qaunt","median","T.Qaunt","max","mean","sd"), las=1, cex.axis = 0.7,
  #         xlab= c(paste("Sammary statistics for the log2norm Table",parameters$recountID)))
  # silence <- dev.off()
  #
  # file.prefix <- paste("diff_bet_stat.log2norm_", parameters$recountID, sep = "")
  # boxplot.file <- file.path(parameters$dir$NormalizationImpact, paste(file.prefix, "boxplot.pdf", sep = "_"))
  # pdf(file = boxplot.file)
  # boxplot(cbind(stat.log2norm$fQuan - stat.log2norm$min, stat.log2norm$median - stat.log2norm$min, stat.log2norm$tQaun - stat.log2norm$min, stat.log2norm$max - stat.log2norm$tQaun) , ylim= c(0,100),
  #         names= c("fQaun-min", "median-min", "tQuan-min", "max-tQuan"), cex.axis = 0.7,
  #         xlab= c(paste("difference between some statistics\n","in log2norm table" ,parameters$recountID)))
  # silence <- dev.off()

  ##### cheking the No. libsum for each classes in the log2norm Table  ##################

  # data.frame(libsum=apply(loaded$dataTable, 1, sum), class=loaded$classes)

  # x <- data.frame(q3 = apply(log2norm$Counts, 1, quantile, q=0.75), sum = apply(log2norm$Counts, 1, sum), class=loaded$classes)
  #
  # x <- data.frame(q3 = apply(rawCounts1, 1, quantile, q=0.75), sum = apply(rawCounts1, 1, sum), class=loaded$classes)

  # x <- data.frame(libsum=apply(loaded$log2norm$counts, 2, sum), class=loaded$log2norm$classLabels)
  #
  # head(x)

  #
  #
  # file.prefix <- paste("libsum.log2norm_",parameters$recountID,sep = "")
  # boxplot.file <- file.path(parameters$dir$NormalizationImpact, paste(file.prefix,"boxplot.pdf", sep = "_"))
  # pdf(file = boxplot.file)
  # save.margins <- par("mar")
  # par(mar = c(14,5,5,1))
  # boxplot(libsum ~ class , data=x, horizontal=TRUE, las=1, cex.axis=0.7,
  #         xlab= c(paste("Libsum statistics for the log2norm Table\n","distribution for each class in",parameters$recountID)))
  # par(mar=save.margins)
  # silence <- dev.off()


  #############################################################################################
  ###### Draw some plots to showing effects of some treatments(Norm and log2) on the Count Table ####
  #############################################################################################
  message.with.time("\tDrawing plots describing count table statistics and distribution")

  #### Compute a trimmed mean: suppress the 5% top and bottom values ####
  if (parameters$compute) {
    message.with.time("Computing trimmed mean of normalized counts")
    x <- unlist(self[["dataTable"]])
    q0.05 <- quantile(x = x, probs = 0.05, na.rm=TRUE)
    q0.95 <- quantile(x = x, probs = 0.95, na.rm=TRUE)
    trimmed <- (x[x > q0.05 & x < q0.95])
    suppressed.proportion <- 1 - length(trimmed)/length(x)
  } else {
    message.with.time("Skipping a trimmed mean")
  }


  ## histogram of normalized counts. Zoom on the representative part of the histogram
  file.prefix <- file.path(plot.dir, paste( dataType ,"_",self[["ID"]], "_","histogram" ,extension, sep = ""))

  pdf(file=file.prefix,
      width = 8, height = 8)
  hist(unlist(self[["dataTable"]]), breaks=10, las=1,
       xlab="normalized counts", ylab="Occurrences", col="grey",
       main=paste(parameters$recountID, "Normalised count distrib"))
  abline(v=mean(unlist(self[["dataTable"]])), lwd=1, col="darkgreen") # mark the mean
  abline(v=median(unlist(self[["dataTable"]])), lwd=1, col="blue") # mark the median
  abline(v=mean(unlist(trimmed)), lwd=1, col="purple") # mark the trimmed mean
  legend("topright",lwd=2,
         legend=c("mean", "median", "trimmed mean"),
         col=c("darkgreen", "blue", "purple"))
  silence <- dev.off()

  #
  # #### Compute a trimmed mean: suppress the 5% top and bottom values ####
  # if (parameters$compute) {
  #   message.with.time("Computing trimmed mean of log2normalized counts")
  #   x <- unlist(self[["dataTable"]])
  #   q0.05 <- quantile(x = x, probs = 0.05, na.rm=TRUE)
  #   q0.95 <- quantile(x = x, probs = 0.95, na.rm=TRUE)
  #   log2.trimmed <- (x[x > q0.05 & x < q0.95])
  #   suppressed.proportion <- 1 - length(log2.trimmed)/length(x)
  # } else {
  #   message.with.time("Skipping a trimmed mean")
  # }
  #
  # file.prefix <- file.path(plot.dir, dataType,paste(self[["ID"]], extension , sep = ""))
  #
  # pdf(file=file.prefix,
  #     width = 8, height = 8)
  # hist(unlist(loaded$log2norm$counts), breaks=100, las=1,
  #      xlab="log2 Normalized counts", ylab="Occurrences", col="grey",
  #      main=paste(parameters$recountID, " log 2 Normalised count distribution"))
  # abline(v=mean(unlist(loaded$log2norm$counts)), lwd=1, col="darkgreen") # mark the mean
  # abline(v=median(unlist(loaded$log2norm$counts)), lwd=1, col="blue") # mark the median
  # abline(v=mean(unlist(log2.trimmed)), lwd=1, col="purple") # mark the trimmed mean
  # legend("topright",lwd=2,
  #        legend=c("mean", "median", "trimmed mean"),
  #        col=c("darkgreen", "blue", "purple"))
  # silence <- dev.off()
  #

  ######  interpretation the difference between the mean and the third quartile raw count #####
  file.prefix <- file.path(plot.dir, paste(dataType, "_",self[["ID"]],"", extension, sep = ""))
  pdf(file= file.prefix,
      width = 8, height = 8)
  plot(x=apply(self[["dataTable"]], 1, mean),
       y=signif(digits=3, apply(self[["dataTable"]], 1, quantile, 0.75)),
       main="raw counts: Percentile 75 versus mean",
       xlab="Mean counts per sample",
       ylab="Percentile 75",
       col=self$sampleColors,
       las=1,
       panel.first=grid())
  legend("topleft", legend=self$classNames, col=self$sampleColors, pch=1, cex=0.8)
  silence <- dev.off()

  ######  interpretation the difference between the mean and the third quartile log2counts #####
  file.prefix <- file.path(plot.dir, paste(dataType,"_",self[["ID"]], "_","mean_3Q",extension, sep = ""))
  pdf(file= file.prefix,
      width = 8, height = 8)
  plot(x=apply(self[["dataTable"]], 1, mean),
       y=signif(digits=3, apply(self[["dataTable"]], 1, quantile, 0.75)),
       main="log2norm counts: Percentile 75 versus mean",
       xlab="Mean counts per sample",
       ylab="Percentile 75",
       col=self$sampleColors,
       las=1,
       panel.first=grid())
  legend("topleft", legend = unique(names(self$sampleColors)), col=self$sampleColors, pch=1, cex=0.8)
  silence <- dev.off()

  message.with.time(" finishing from drawing plots describing count table statistics and distribution")

}
