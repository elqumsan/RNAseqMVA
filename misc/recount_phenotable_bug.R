#### Genes as features: the phenotable extraction works fine ####

## Download recount data if not already there
data.dir <- "~/RNAseqMVA_workspace/data/SRP056295/"
dir.create(data.dir, recursive = TRUE, showWarnings = FALSE)
setwd(data.dir)
if (!file.exists("rse_gene.Rdata")) {

}

## Load a recount memory image
load("rse_gene.Rdata")
gene_counts <- scale_counts(rse_gene, by  = "mapped_reads")
gene_pheno <- colData(gene_counts)
gene_geochar <-   recount::geo_characteristics(gene_pheno)
names(gene_geochar)
head(gene_geochar)
table(gene_geochar$tissue)

#### Transcripts (tx) as features: the phenotable extraction fails ####

## Load a recount memory image
load("rse_tx.RData")
tx_counts <- scale_counts(rse_tx, by  = "mapped_reads")
tx_pheno <- colData(tx_counts)
tx_geochar <-   recount::geo_characteristics(tx_pheno)

names(gene_geochar)
## Problem: there is a single column named "c..tissue" instead of 2 columns named "tissue" and "cell type"

head(tx_geochar)
## This apparently comes from a parsing problem: the column name (cell type)
## comes in the column content rather than as a column header


