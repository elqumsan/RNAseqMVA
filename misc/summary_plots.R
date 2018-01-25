
#### Select a combination of recountIDs, classifiers and suffixes ####
tests <- expand.grid(
  classifier = c("knn", "rf", "svm"),
  recountID = recountIDs <- c("SRP042620",
                              "SRP056295",
                              "SRP061240",
                              "SRP062966"),
  data.type = "log2norm",
  suffix = c("allvars.tsv")
)


## Tester
tests$files <- file.path("results",
                   tests$recountID,
                   tests$classifier,
                   "tables",
                   paste(sep="_", tests$classifier,
                         #                                "_", tests$recountID,
                         tests$data.type,
                         tests$suffix))
tests$perm.files <- sub(pattern = ".tsv", replacement = "_permLabels.tsv", x = tests$files)

#### Select experiments for which both real data and permuted labels files exist ####
tests$both.files.ok <- file.exists(tests$files) & file.exists(tests$perm.files)
# View(tests)

#### Read error tables ####
error.summary <- SummarizeErrorTable(tests$files, stopIfMissing = FALSE)
perm.error.summary <- SummarizeErrorTable(tests$perm.files, stopIfMissing = FALSE)

#### Plot comparing  error rates with real and permuted labels ####
#names(tests)
colors.per.classifier <- c("knn" = "blue",
                           "svm" = "brown",
                           "rf" = "darkgreen")
tests$color <- colors.per.classifier[as.vector(tests$classifier)]

pch.per.recountID <- c("SRP042620" = 1,
                       "SRP056295" = 2,
                       "SRP061240" = 3,
                       "SRP062966" = 4)
tests$pch <- pch.per.recountID[as.vector(tests$recountID)]

tests[, c("classifier", "color")]

plot(error.summary$testing.error.rate,
     perm.error.summary$testing.error.rate,
     xlim = c(0,1), ylim = c(0,1),
     xlab = "Error rate",
     ylab = "Error rate, permuted lables",
     col = tests$color,
     pch = tests$pch,
     panel.first = c(grid(),
                     abline(a=0,b=1, col="grey")))
legend("bottomright",
       legend=names(pch.per.recountID),
       pch=pch.per.recountID)
legend("topright",
       legend=names(colors.per.classifier),
       col=colors.per.classifier,
       pch=1)

