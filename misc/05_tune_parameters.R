#### Use e1071::tune() wrapper to find optimal parameters for the different classifiers ####

## Note: e1071::tune() functions are versy convenient to analyse user-specified combinations of
## parameters but they run all the tests sequentially, which may take a huge time for big datasets
## (as the ones we analyse) and if the number of parameter values is high.
##
## Moreover, the tuning tests each parameter value only once, but we noticed that for some parameters
## there are strong fluctuations between successive runs (e.g. the ntrees and mtry parameters of randomForest).
## For this reason, we developed a specific procedure to run evaluate the impact of parameters with a
## user-specified number of iterations (which enables us to evaluate the fluctuatins of performances
## between independent runs).
##
## We however keep this code for the sake of comparison.

## Run the analysis

