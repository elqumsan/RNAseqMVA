
################################################################
## Run differential analysis with DESeq2 and edgeR to define variable order


## Run differential analysis with edgeR to define variable order
message.with.time("Running edgeR to define variable ordering")
DEG.edgeR  <- DEGordering(loaded$countTable, loaded$classes, method = "edgeR")
message.with.time("Done edgeR DEG analysis")
# sorted.log2.transformed.edgeR <- loaded$countTable[, DEG.DESeq2$geneOrder]

## Run differential analysis with DESeq2 to define variable order
message.with.time("Running DESeq2 to define variable ordering")
DEG.DESeq2 <- DEGordering(loaded$countTable, loaded$classes, method = "DESeq2")
message.with.time("Done DESeq2 DEG analysis")
# sorted.log2.transformed.DESeq2 <- loaded$countTable[, DEG.edgeR$geneOrder]

## Order genes at random
message.with.time("Generating a fake DEG result by random gene ordering")
DEG.randomized <- list()
DEG.randomized$geneOrder <- sample(DEG.edgeR$geneOrder)

message.with.time("Finishing from the DEG Computation Process")
