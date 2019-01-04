#### Bug with the pheno table of transcripts ####
# There is a bug with the pheno tables of the rse_tx objects.
# For several experiments, the "characteristics" column of the DataFrame returned by colData(rse_tx) contains strangely placed quotes which perturb the parsing.
# I paste below a minimal code that reproduces the bug.
# The same bug occurs with several recount IDs, but not all.


#### Gene-wise counts (this first part works fine) ####

## Download data in rse-gene format
recountID <- "SRP056295"
gene_url <- download_study(project = recountID, type = "rse-gene", download = TRUE)
print(gene_url)

## Load the rse_gene object in memory
load(file.path(recountID, 'rse_gene.Rdata'))

## Extract GEO characteristics from the rse_gene object
gene_geochar <- recount::geo_characteristics(colData(rse_gene))
head(gene_geochar)
table(gene_geochar)

#### Transcript-wise counts #####

## Download the rse-tx object
tx_url <- download_study(project = recountID, type = "rse-tx", download = TRUE)
print(tx_url)

## Inconsistency: the following line fails on Linux systems because the extension
## is RData for transcripts, whereas it is Rdata for genes.
## It works on Mac OS X because the system is flexible with file upper/lower cases.
load(file.path(recountID, 'rse_tx.Rdata'))

## This works on Linux as well as Mac OS X
load(file.path(recountID, 'rse_tx.RData'))

## Extract GEO characteristics from the rse_gene object
tx_geochar <- recount::geo_characteristics(colData(rse_tx))
head(tx_geochar)
table(tx_geochar)

## The bug apparently comes from the pheno table
head(colData(rse_tx)$characteristics)


## Fruitless attempt to debug
pheno <- colData(rse_tx)
head(pheno$characteristics)
pheno$characteristics <-
  gsub(
    pattern = '\\"', perl = TRUE,
    replacement = "'",
    x = pheno$characteristics)
head(pheno$characteristics)

tx_geochar2 <- recount::geo_characteristics(pheno)
head(tx_geochar2)
table(tx_geochar2)

