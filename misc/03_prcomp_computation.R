
#### Compute principal components for normalized log2 counts ####
message.with.time("Pre-processing by Principal Component analysis (PCA)")
log2norm.prcomp.centred <- prcomp( na.omit(log2norm), center = TRUE, scale. = FALSE)
log2norm.prcomp.centred.scaled <- prcomp(na.omit(log2norm), center = TRUE, scale. = TRUE)

## Indicate that this principal components analysis for log2 count has finished running
message.with.time("finished running Principal Component analysis (PCA) for normalized log2 counts")

##### Exhibiting the geo charactiristics for the current project #####
message.with.time("Exhibit the geo charactiristics for such experiment: ", parameters$recountID, " in order to know the class lable
                  for such experiment")
head( geo.characteristics)
geo.characteristics.file <- file.path("~/RNAseqMVA_workspace", "data", parameters$recountID, "geo.characteristics.tsv")
write.table( geo.characteristics, file = geo.characteristics.file, quote = FALSE,
             row.names = FALSE, sep = "\t" )
