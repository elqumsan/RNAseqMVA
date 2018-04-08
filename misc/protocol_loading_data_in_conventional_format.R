#' @title a minimal protocol to providing the minimal instructions to load the data in a conventional format (tab-delimited, everything in the form of data.frames)
#' @author  Mustafa ABUELQUMSAN and Jacques van Helden
#' @description for the safe of facilitate the data exchanges with statisticians and biloigest we expoert all the RNAseq data in a conventional formate (tab-delimited, everything in the form of data.frames)

#' @param
#' @return
#' @import
#' @export


# How to laod the raw count table in conventional format as data.frame
rawCount.dataframe <-as.data.frame( read.delim(paste( parameters$dir$tsv,"/rawCounts_", parameters$recountID,".tsv", sep = "")))
dim(rawCount.dataframe)

# How to load the Phento table in conventional format as data.frame
Pheno.dataframe <- as.data.frame(read.delim(paste(parameters$dir$tsv, "/pheno_", parameters$recountID,".tsv", sep = "")))

dim(Pheno.dataframe)
View(Pheno.dataframe)

# How to load Normalized count Table
normCounts.dataframe <- as.data.frame(read.delim(paste(parameters$dir$tsv,"/NormCounts_",parameters$recountID, ".tsv", sep = "")))
dim(normCounts.dataframe)
