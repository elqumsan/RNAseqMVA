
#### Compute principal components for normalized log2 counts ####
message.with.time("Pre-processing by Principal Component analysis (PCA)")
log2norm.prcomp.centred <- prcomp( na.omit(log2norm), center = TRUE, scale. = FALSE)
log2norm.prcomp.centred.scaled <- prcomp(na.omit(log2norm), center = TRUE, scale. = TRUE)

## Indicate that this principal components analysis for log2 count has finished running
message.with.time("finished running Principal Component analysis (PCA) for normalized log2 counts")


