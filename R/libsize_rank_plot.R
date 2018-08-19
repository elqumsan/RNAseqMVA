#' @title Rank plot of RNA-seq library sizes
#' @author Jacques van Helden
#' @description Plot the sizes (in million reads) of RNA-seq librairies ranked by decreasing size.
#' @param count.table a table of counts per feature, with one row per feature (e.g. gene, transcript) and one column per library (e.g. run, sample).
#' @param plot.file if ot null, save th plot in pdf format in the specified file
#' @param plot.width=7 width of the pdf file (only works when plot.file is not null)
#' @param plot.height=5 height of the pdf file (only works when plot.file is not null)
#'
#' @export
LibsizeRankPlot <- function(count.table,
                            plot.file = NULL,
                            plot.width = 7,
                            plot.heigh = 5
                            ) {
  message("\tPlotting ranked libsizes")
  message("\t\t", plot.file)
  if (!is.null(plot.file)) {
    pdf(file = plot.file, width = plot.width, height = plot.heigh)
  }
  filtered.libsizes <- apply(count.table, 2, sum) / 1e+6
  plot(sort(filtered.libsizes, decreasing = TRUE),
       type = "h", col = "grey", las = 1,
       xlab = "Sample rank",
       ylab = "sample size (Mreads)",
       ylim = c(0, ymax),
       main = paste(sep = "", recountID, " \n", "library sizes"))
  ymin <- min(filtered.libsizes)
  ymax <- max(filtered.libsizes)
  abline(h = ymin, col = "orange", lwd = 2)
  abline(h = ymax, col = "#22BB44", lwd = 2)
  arrows(length = 0.1, angle = 15, lwd = 2, col = "red", code = 3,
         x0 = length(filtered.libsizes),
         x1 = length(filtered.libsizes),
         y0 = ymin,
         y1 = ymax)

  text(x = length(filtered.libsizes), y = (ymax + ymin) / 2,
       labels = paste(round(ymax / ymin * 100), "%"),
       pos = 2)
  if (!is.null(plot.file)) {
    silence <- dev.off(); rm(silence)
  }

}
