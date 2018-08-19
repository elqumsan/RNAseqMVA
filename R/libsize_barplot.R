#' @title Draw a barplot with the number of million reads per sample
#' @author Jacques van Helden
#' @param counts a count table with 1 row per feature and 1 column per sample
#' @param sample.labels=colnames(counts) sample labels (displayed as y legend)
#' @param sample.colors=NULL sample-specific colors (e.g. reflecting conditions)
#' @param main="Library sizes" Main title for the plot
#' @param xlab="libsum (Million reads per sample)" X axis label
#' @param plot.file=NULL Path of a file to store the figure (pdf format).
#' @param fig.format="pdf" File format (only used if plto.file is defined)
#' @param cex.names=1 fondt size for sample labels
#' @param lmargin=cex.names*max(nchar(sample.labels))/2.5+5 Adapt margin to the label sizes
#' @param width=7 width of the plot file
#' @param height=length(sample.labels)/3+2 height of the plot file (if not specified, will be adapted to the number of samples)
#' @param ... additional parameters are passed to barplot
#' @export
LibsizeBarplot <- function(counts,
                           sample.labels=colnames(counts),
                           sample.colors=NULL,
                           main = "Library sizes",
                           xlab = "libsum (Million reads per sample)",
                           # plot.file = NULL, fig.format = "pdf", 
                           cex.names = 1,
                           lmargin = cex.names * max(nchar(sample.labels))/2.5 + 5,
                           # width = 7,
                           # height = length(sample.ids)/3 + 2,
                           ...) {

  ## Make sure that sample colors and sample labels are provided as a vector rather than a factor or a list
  sample.labels <- as.vector(unlist(sample.labels))
  sample.colors <- as.vector(unlist(sample.colors))
  
  ## Store original graphical parameters
  par.ori <- par(no.readonly = TRUE)

  ## Compute library sizes
  libsizes <- round(x = apply(counts, 2, sum, na.rm = TRUE) / 1e6, digits = 2)

  ## Sample colors
  if (is.null(sample.colors)) {
    #sample.colors <- colorRampPalette(c('blue','red'))(100)[as.numeric(cut(libsizes, breaks = 100))]
    sample.colors <- rep(c("#888888", "#DDDDDD"), length.out = length(libsizes))
  }

  ## Sample-wise library sizes
  # if (!is.null(plot.file)) {
  #   message("Generating plot", plot.file)
  #   file.prefix <- sub(x = plot.file,
  #                      pattern = paste(sep = "", ".", img.format), replacement = "")
  #   
  #   OpenPlotDevice(file.prefix = file.prefix, fig.format = fir.format, width = width, height = height)
  # }

  par(mar = c(5, lmargin, 4, 1)) ## adapt axes
  bplt <- barplot(libsizes,
                  names.arg = sample.labels,
                  main = main,
                  horiz = TRUE, las = 1,
                  xlab = xlab,
                  col = sample.colors, ...)
  grid(col = "white", lty = "solid", ny = 0)
  text(x = pmax(libsizes, 3),
       labels = libsizes, y = bplt, pos = 2, font = 2)
  # if (!is.null(plot.file)) {
  #   silence <- dev.off(); rm(silence)
  # }
  par(par.ori)
#  return(bplt)
}

