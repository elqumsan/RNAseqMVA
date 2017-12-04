#' @title Select Recount projects based on the presence of user-specified fields in their characteristics.
#' @description Recount metadata contains a field named "characteristics", which is a list of fields with values.
#' It is useful to run a query on the field names (not values) in order to select projects
#' related to a given experimental setting (e.g. "time series"), or to a particular field (e.g. "disease type").
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param recount.metadata a data.frame containing the metadata, as returned by recount::all_metadata() function.
#' @param query.charact a string indicating the query to be matched against field names in the characteristics tables.
#' @param run.kable=FALSE if TRUE, print out the summary of the matches in kable format for markdown report.
#' @param run.barplot=FALSE if TRUE, draw a barplot with the number of runs per selected project
#' @examples
#' ## Load recount metadata
#' library(recount)
#' recount.metadata <- all_metadata(subset='sra')
#'
#' ## Select all the projects having "disease status" as a characteristics field
#' selectProjectsByCharacteristics(recount.metadata,
#'    query.charact = "disease status",
#'    run.barplot = TRUE)
#' @import recount
#' @export
selectProjectsByCharacteristics <- function(recount.metadata,
                                            query.charact,
                                            run.kable = FALSE,
                                            run.barplot = FALSE) {

  result <- list()
  result$query.charact <- query.charact

  ## Concatenate characteristics (list of fields) in a single string per run
  characteristics.string <- unlist(lapply(recount.metadata$characteristics, paste, collapse="; "))

  ## Back-conversion of "NA" strings (obtained from NA values) to NA values
  characteristics.string[characteristics.string == "NA"] <- NA

  ## Count the number of NA or non-NA values for characteristics
  charact.defined <- as.data.frame.table(table(!is.na(characteristics.string)))
  names(charact.defined) <- c("Characteristics defined", "Number of runs")

  ## Display a table with the number of runs having defined characteristics attributes
  if (run.kable) {
    kable(charact.defined, caption = "Number of runs with defined or undefined characteristics in Recount2 metadata. ")
    #table(is.na(unlist(recount.metadata$characteristics)))
  }

  ## Cast metadata to a data frame with a single row per run
  ## Use the concatenated chacteristics to fill the column "characteristics"
  recount.metadata.frame <- as.data.frame(recount.metadata)
  recount.metadata.frame$characteristics <- characteristics.string
  # table(is.na(recount.metadata.frame$characteristics))

  ## Return all projects havin a field named "disease status"
  matching.runs <- grep(recount.metadata.frame$characteristics, pattern = query.charact)
  result$matching.runs <- matching.runs
  if (length(matching.runs) == 0) {
    message("Not a single match in recount characteristics fields for the query ", query.charact)
    return(result)
  }

  ## Get the project ID associated to each run
  matching.project.id.per.run <- recount.metadata.frame[matching.runs, "project"]
  result$matching.project.id.per.run <- matching.project.id.per.run

  ## Get a unique list of matching project IDs
  matching.project.ids <- unique(sort(matching.project.id.per.run))
  result$matching.project.ids <- matching.project.ids

  ## Count the number of runs per matching project
  matching.projects.Nruns <- as.data.frame.table(sort(table(matching.project.id.per.run), decreasing=TRUE), row.names = 1)
  names(matching.projects.Nruns) <- "Nruns"
  result$matching.projects.Nruns <- matching.projects.Nruns

  if (run.kable) {
    kable(matching.projects.Nruns, caption=paste("Recount projects having", query.charact, "as field name for the characteristics. "))
  }

  if (run.barplot) {
    par.ori <- par(no.readonly = TRUE)
    par(mar=c(4,6,5,1))
    barplot(as.vector(as.matrix(matching.projects.Nruns)),
            names.arg = rownames(matching.projects.Nruns),
            main=paste("Query field: query.charact"),
            horiz = TRUE, las = 1, cex.names = 0.7, cex.axis = 0.8,
            xlab="Runs per project", col="#88FFDD")
    legend("topright",
           c(paste(length(matching.project.ids), "projects"),
             paste(length(matching.project.id.per.run), "runs")))
    par(par.ori)
  }
  return(result)
}
