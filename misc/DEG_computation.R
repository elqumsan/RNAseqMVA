
################################################################
#### Run differential analysis with DESeq2 and edgeR to define variable order ####


#### Run differential analysis with edgeR to define variable order ####
message.with.time("Running edgeR to define variable ordering")
DEG.edgeR  <- DEGordering(loaded$countTable, loaded$classes, method = "edgeR")
message.with.time("Done edgeR DEG analysis")
# sorted.log2.transformed.edgeR <- loaded$countTable[, DEG.DESeq2$geneOrder]

#### Run differential analysis with DESeq2 to define variable order ####
message.with.time("Running DESeq2 to define variable ordering")
DEG.DESeq2 <- DEGordering(loaded$countTable, loaded$classes, method = "DESeq2")
message.with.time("Done DESeq2 DEG analysis")
# sorted.log2.transformed.DESeq2 <- loaded$countTable[, DEG.edgeR$geneOrder]

#### Order genes at random ####
message.with.time("Generating a fake DEG result by random gene ordering")
DEG.randomized <- list()
## such is randomizing the ordered count table by DEG edgeR; BUT be ware we should keep the neme of this variable as
## DEG.randomized$orderedCountTable to we can correctly use it in the script (DEG_impact_of_top_gene_number.R) where there
## is method to automatically get the name of object under the "$orderedCountTable"; so that you should keep this part
# as it that to avoiding the vey critical error in computation.
DEG.randomized <- DEG.edgeR
DEG.randomized$geneOrder <-sample( DEG.randomized$geneOrder)
DEG.randomized$orderedCountTable <- DEG.edgeR$orderedCountTable[ ,sample(DEG.edgeR$geneOrder)]

message.with.time("Finishing from the DEG Computation Process")
