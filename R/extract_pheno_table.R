## Mustafa will adda  proper roxygen2 description


ExtractPhenoTable <- function(rse,
                              verbose = FALSE) {

  result <- list()

  ## Table with information about the columns of the RangedSeummaryExperiment.
  if (verbose) {
    message("Building pheno table")
  }
  runPheno <- colData(rse) ## phenotype per run
  # pheno <- runPheno ## A TRICK
  # View(runPheno)
  # names(runPheno)


  ## Extract the conditions from the "characteristics" column of the coldata.
  ## This is a bit tricky: we have to parse a string describing several attributes.
  geochar <- lapply(split(
    runPheno,
    seq(from=1, to=nrow(runPheno))),
    geo_characteristics)
  # View(geochar)
  # names(geochar)
  # head(geochar)

  geochar <- do.call(rbind, lapply(geochar, function(x) {
    if('cells' %in% colnames(x)) {
      colnames(x)[colnames(x) == 'cells'] <- 'cell.line'
      return(x)
    } else {
      return(x)
    }
  }))

  # View(geochar)
  # head(geochar)

  ## Build a pheno table with selected columns from coldata + the geodata we just extracted
  runPhenoTable <- cbind(
    runPheno[, grep(pattern="(characteristics|sharq)", x=names(runPheno), invert=TRUE)],
    geochar)
  # View(phenoTable)
  # class(phenoTable)

  ## Extract a phenoTable with selected fields from the runPheno object
  runPhenoTable2 <- data.frame(
    project = runPheno$project,
    sample = runPheno$sample,
    experiment = runPheno$experiment,
    run = runPheno$run,
    geo_accession = runPheno$geo_accession,
    characteristics = runPheno$characteristics@unlistData
  )

  ## Missing: parse sub-fields from the "characteristics" field (somehow tricky)

  rownames(runPhenoTable2) <- runPhenoTable2$run
  # View(runPhenoTable2)
  # class(runPhenoTable2)
  result$runPhheno <- runPheno
  result$geochar <- geochar
  result$runPhenoTable <- runPhenoTable
  result$runPhenoTable2 <- runPhenoTable2
  return(result)
}
