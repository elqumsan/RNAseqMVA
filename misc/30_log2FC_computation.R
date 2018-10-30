## Computing log2FC
## For didactic purposes


recountID <- "SRP048759"
studyCases <- loadCounts(recountID = recountID, mergeRuns = T, classColumn = "tissue")

## We choose  TMM as normalization method because it gave good results in Dillies et al. (2013)
log2normCounts <- studyCases[[recountID]]$datasetsForTest$TMM_log2

names(log2normCounts)

counts <- log2normCounts$dataTable
classLabels <- log2normCounts$classLabels
classes <- unique(classLabels)

## Let us assume that we want to compare reference and test conditions
ref.condition <- classes[1]
test.condition <- classes[2]

ref.counts <- counts[, classLabels == ref.condition]
dim(ref.counts)
ref.mean.per.gene <- apply(ref.counts, 1, mean)

test.counts <- counts[, classLabels == test.condition]
dim(test.counts)
test.mean.per.gene <- apply(test.counts, 1, mean)

## XY plot of mean counts per gene for test condition versus ref condition
plot(ref.mean.per.gene, test.mean.per.gene,
     main = "Mean log2 normalized counts per gene",
     pch = 1,
     col = densCols(x = ref.mean.per.gene, y = test.mean.per.gene),
     xlab = ref.condition,
     ylab = test.condition)
grid()
abline(a = 0, b = 1, col = "black")

## MA plot :
## - abcissa represents the mean expression value (across all samples of ref + test)
## - ordinate represents the log2FC

## Since the counts are already log2-transformed by our normalization procedure
## we compute the difference between means of log2-coutns, which corresponds to the log2 of the ration between counts/
M <-  test.mean.per.gene - ref.mean.per.gene
A <- apply(cbind(test.counts, ref.counts), 1, mean) ## Average of log2-transformed counts

plot(x = A, y = M,
     main = "(approx) MA plot",
     pch = 1,
     col = densCols(x = A, y = M),
     xlab = ref.condition,
     ylab = test.condition)

grid()
abline(h = 0, col = "black")
