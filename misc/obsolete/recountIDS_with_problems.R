#### RecountIDs with problems (TO INVESTIGATE / DEBUG LATER) ####
recountIDs.with.problems <- c("SRP003611" = "the number of samples drops after run merging",
                              "SRP039694" = '	Building class-specific attributes for DataTableWithClasses	SRP039694
Error in `colnames<-`(`*tmp*`, value = c("Class", "nbSamples")) :
  "names" attribute [2] must be the same length as the vector [1]',
                              "SRP008976" = "",
                              "SRP006574" = "	Building pheno table\nError in rbind(deparse.level, ...) : numbers of columns of arguments do not match",
                              "SRP042161" = "class column not specified",
                              "SRP041736" = "We cannot analyse this dataset because the pheno table does not contain any info about the sample classes
class column not specified (NA)
Instead of complaning because of this, the program crashes with the following message:

TO DO: detect such cases and stop with explicit message before any analysis.

2018-04-07_093304		Creating object of class DataTableWithClasses
	Building attributes for object of class DataTableWithClasses
                              Building class-specific attributes for DataTableWithClasses	SRP041736
                              Error: subscript contains invalid names
                              ")
