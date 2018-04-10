#' @title Plot two components of prcomp transform
#' @author Jacques van Helden
#' @description Plot the first and second components of a prcomp() result in a lighter way than the default biplot.
#' @param self an object of class countTableWithClasses containing a field prcomp (produced by prcompt function)
#' @param main main title is built by default by can be replaced
#' @param col = self$sampleColors dots are colored according to sample class membership.
#' @param pcs=c(1,2) the two indices of the principal components to be displayed
#' @param ... additional parameters are passed to plot()
#' @export
#'
plot2PCs <- function(self,
                     main = paste(sep="", self$ID, "; PCs of log2norm counts"),
                     col = self$sampleColors,
                     pcs = c(1,2),
                     ...) {



  ## Check conformity of the arguments
  if (!is(object = self, class2 = "countTableWithClasses")) {
    stop("plot2PCs(). Invalid self argument for plotPC12(): must be an object of class countTableWithClasses. ")
  }
  if (is.null(self$prcom)) {
    stop("plot2PCs(). This object does not contain a prcomp attribute)")
  }
  if (!is(object = self$prcomp, class2 = "prcomp")) {
    stop("plot2PCs(). prcomp attribute of self must belong to prcomp class.")
  }

  pcs <- as.integer(pcs)
  if ((length(pcs) != 2) ||
      (!is.integer(pcs)) ||
      (min(pcs) < 1) ||
      (max(pcs) > ncol(self$prcomp$x))) {
    stop("plot2PCs(). Argument pcs should contain 2 integer values comprised between 1 and the number of components. ")
  }

  ## Get the two first components for the plot
  PC1 <- self$prcomp$x[,pcs[1]]
  PC2 <- self$prcomp$x[,pcs[2]]

  ## Compute required space for the legend according to the number of classes
  ylim <- c(min(PC2), max(PC2) + self$nbClasses*0.05*(diff(range(PC2))))

  class.pch <- 1:self$nbClasses
  names(class.pch) <- self$classNames
  sample.pch <- class.pch[self$classLabels]

  ## Plot PC1 vs PC2
  plot(x = PC1,
       y = PC2,
       main = main,
       col = col,
       xlab=paste("PC", pcs[1]),
       ylab=paste("PC", pcs[2]),
       ylim = ylim, las=1,
       pch = sample.pch,
       panel.first = c(
         grid(col="#DDDDDD", lty = "solid"),
         abline(v=0, col="black"),
         abline(h=0, col="black")),
       ...
  )
  legend("topleft",
         legend = studyCases[[recountID]]$log2normPCs$classNames,
         col = studyCases[[recountID]]$log2normPCs$classColors,
         border="black", bty='o', bg='white',
         pch = class.pch, cex = 0.8
  )


}
