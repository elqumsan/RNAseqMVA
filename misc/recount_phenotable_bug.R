## This scripts aims at testing the loading of phenoTable with the recount package,  for genes and transcripts respectively

#### Load recount library and check version ####
library(recount)
# message(version$version.string)
message("Recount version: ", packageDescription("recount")$Version)

#### Re-install recount package if old version ####
if (packageDescription("recount")$Version < "1.10.13") {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("recount")
}

#### Parameters ####
recountID <- "SRP056295" ## an example of problematic dataset
data.dir <- file.path("~/RNAseqMVA_workspace/data/", recountID)

#### Download recount data if not already there ####
dir.create(data.dir, recursive = TRUE, showWarnings = FALSE)
setwd(data.dir)
list.files()
if (!file.exists("rse_gene.Rdata")) {
  gene_url <- download_study(project = recountID, type = "rse-gene", download = TRUE)
  message("\tdownloaded gene Rdata from ", gene_url)
}
if (!file.exists("rse_tx.RData")) {
  tx_url <- download_study(project = recountID, type = "rse-tx", download = TRUE)
  message("\tdownloaded transcript Rdata from ", tx_url)
}

#### Genes as features: the phenotable extraction works fine ####
message("Loading gene data")
load("rse_gene.Rdata")
gene_counts <- scale_counts(rse_gene, by  = "mapped_reads")
gene_pheno <- colData(gene_counts)
gene_geochar <-   recount::geo_characteristics(gene_pheno)
names(gene_geochar)
print(head(gene_geochar))
print(table(gene_geochar$tissue))

#### Transcripts (tx) as features: the phenotable extraction fails ####

## Load a recount memory image
message("Loading transcript data")
load("rse_tx.RData")
tx_counts <- scale_counts(rse_tx, by  = "mapped_reads")
tx_pheno <- colData(tx_counts)
tx_geochar <-   recount::geo_characteristics(tx_pheno)

names(gene_geochar)
print(head(tx_geochar))
print(table(tx_geochar$tissue))


