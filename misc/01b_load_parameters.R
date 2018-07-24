require("RNAseqMVA")

configFile <- "~/RNAseqMVA/misc/00_project_parameters.yml"
message.with.time("Loading parameters from YAM file ", configFile)

## Read global and study case-specific parameters from a yaml-formatted file.
## These parameters will then be used to overwrite the default parameters.
project.parameters <- yaml.load_file(configFile)

selectedRecountIDs <- project.parameters$global[["selected_recount_ids"]]

## Get all recount IDs
recountIDs <- grep(pattern = "^SRP", x = names(project.parameters), value = TRUE)
message("\tYAML config file contains ", length(recountIDs)," recount IDs.")

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

## Optional: select a subset of the recountIDs
# selectedRecountIDs <- c("SRP057196") ## single-cell of vairous cell types in different brain tissues
#selectedRecountIDs <- c("SRP056295") ## Human leukemia
# selectedRecountIDs <- c("SRP042620")
#selectedRecountIDs <- c("SRP042620", "SRP057196") #, "SRP056295")
#selectedRecountIDs <- setdiff(recountIDs, names(recountIDs.with.problems))



message("\tSelected ", length(selectedRecountIDs)," recount IDs: ", paste(collapse = "; ", selectedRecountIDs))
#selectedRecountIDs <- recountIDs

# recountID <- "SRP042620"   ## Multi-group breast cancer
# recountID <- "SRP057196"    # individual cells into all of the major neuronal, glial, and vascular cell types in the brain

# recountID <- "SRP003611"   # transcriptomes of 347 cells from 10 distinct populations in both of low-coverage (~0.27 million reads per cell) and high-coverage (~5 million reads per cell)
# recountID <- "SRP061240"   # several types from cancer (pancreatic, colorectal, prostat cancer) against Healthy control

#recountID <- "SRP006574"  # NOT working:  # this Project will constitute the classes that have IRanges object and then it will not have the same length for count Table such leads for STOP

# recountID <- "SRP062966"   # this project Boold disease
# recountID <- "SRP066834"   # cerebral cancer type.
# recountID <-  "SRP039694"  # hepatocellular carcinoma
# recountID <- "SRP008976"   # NOT working properly
######################### such these projects have promptly working #####################
# recountID <- "SRP042620"   ## Multi-group breast cancer. working excellent.
# recountID <- "SRP061240"   # several types from cancer (pancreatic, colorectal, prostat cancer) against Healthy control
# recountID <- "SRP062966"   # this project Blood disease. working excellent

# recountID <- "SRP056295"  # types of leukemia there is problem with random forest it say there is missing classes BUT actually there aren't any missing in the y class lable

# recountID <- "SRP003611"   # transcriptomes of 347 cells from 10 distinct populations in both of low-coverage (~0.27 million reads per cell) and high-coverage (~5 million reads per cell)
# recountID <- "SRP006574"   # NOT working:  # this Project will constitute the classes that have IRanges object and then it will not have the same length for count Table such leads for STOP
# recountID <- "SRP066834"   # Not working: cerebral cancer type.
# recountID <-  "SRP039694"  # Not working: hepatocellular carcinoma
# recountID <- "SRP008976"   # NOT working properly
# recountID <- "SRP006575"   # Not working


# message("\tDefault recountID ", recountID)

#### Initialise cluster processing ####
if (!is.null(project.parameters$global$jobs)) {
  ## NOTE: WE SHOULD HAVE A TERMINATING SCRIPT, FOR INSTANCE TO STOP THE CLUSTER
  ## Stop the cluster
  #  stopImplicitCluster()
  project.parameters <- initParallelComputing(project.parameters)
  cl <- project.parameters$global$cl
}


## END OF SCRIPT
#################################################################

