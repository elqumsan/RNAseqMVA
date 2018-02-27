# # Define a directory for PCA figures
# dir.visualisePCs <- file.path(dir.results, "PCA")
# dir.create(dir.visualisePCs, showWarnings = FALSE, recursive = TRUE)

################################################################
## Compute principla components for normalized log2 counts
log2norm.prcomp.centred <- prcomp(log2norm$Counts, center = TRUE, scale. = FALSE)
log2norm.prcomp.uncentred <- prcomp(log2norm$Counts, center = FALSE, scale. = FALSE)
# log2norm.prcomp.centred.scaled <- prcomp(log2normCounts, center = TRUE, scale. = TRUE)
# log2norm.prcomp.scaled <- prcomp(log2normCounts, center = T , scale. = T)

## Explore the principal component results
names(log2norm.prcomp.centred)
head(log2norm.prcomp.centred$center)
head(log2norm.prcomp.centred$scale)
dim(log2norm.prcomp.centred$rotation)
log2norm.prcomp.centred$rotation[1:10, 1:5] ## This corresponds to the loadings, i.e. the weight of each variable (gene) on each component
log2norm.prcomp.centred$x[1:10, 1:5] ## x corresponds to the "scores", or coordinates of each individual (sample) on the components
## The number of dimension has automatically been restricted to the number of samples


## We explored the princomp method but actually it is not appropriate because
## it does not accept data with more variables than individuals, which which
## always ou case.
# log2norm.princomp <- princomp(t(log2normCounts) , cor=FALSE)
# names(log2norm.princomp)
# log2norm.princomp$call
# loadings(log2norm.princomp)
# dim(log2norm.princomp$scores)
# log2norm.princomp$scores[1:10, 1:5]
# log2norm.princomp$loadings[1:10, 1:5]

###### Stadard deviation plot ######
PC.sd.plot <- file.path(dir.visualisePCs, paste("standard_deviation_plot.pdf"))
pdf(file = PC.sd.plot)
barplot(log2norm.prcomp.centred$sdev,
        main= paste0("centred PCs (no scaling) with experiment no.", parameters$recountID),
        las=1, ylab="Standard dev per coponent", xlab="Component number")
silence <- dev.off()

##### Variance plot: compare scaled and unscaled #####
PC.var.plot <- file.path(dir.visualisePCs, paste("compare_scaled_and_unscaled_effect.pdf"))
message("Plotting scaled and unscaled of principal components for the experiment no.", parameters$recountID, " in ",file.path(PC.var.plot))
pdf(file =PC.var.plot )
par(mfrow = c(2,1))
barplot((log2norm.prcomp.centred$sdev^2)[1:50], las=1, ylab="Variance per coponent", xlab="Component number", main="Centered variables")
barplot((log2norm.prcomp.uncentred$sdev^2)[1:50], las=1, ylab="Variance per coponent", xlab="Component number", main="Uncentered variables")
par(mfrow = c(1,1))

silence <- dev.off()

####### Plot pairs of principal components ######
PC.pair.file <- file.path(dir.visualisePCs, paste0("PCA_plot_first_compontent_pairs_",parameters$recountID,".pdf"))
message("Plotting pairs of principal components for the experiment no.", parameters$recountID," file path", PC.pair.file)

pdf(file = PC.pair.file, width=12, height=12)
par(mfrow=c(2,2))
## Plot each sample on the 2 first components
## each dot is one patient
ymin <- PC2.min <- min(log2norm.prcomp.centred$x[,2])
PC2.max <- max(log2norm.prcomp.centred$x[,2])
ymax <- PC2.max + 0.4*(PC2.max - PC2.min)
plot(log2norm.prcomp.centred$x[,c(1,2)],
     col=sampleColors,
     main= paste0("PC plot for the experiment no._",parameters$recountID), las=1, panel.first=grid(),
     ylim=c(ymin, ymax))
legend("topright", legend=names(classColors), col=classColors, cex=1, pch=1)

plot(log2norm.prcomp.centred$x[,c(3,4)],
     col=sampleColors,
     main=paste0("PC plot for the experiment no._", parameters$recountID), las=1, panel.first=grid())

plot(log2norm.prcomp.centred$x[,c(1,3)],
     col=sampleColors,
     main=paste0("PC plot for the experiment no._", parameters$recountID), las=1, panel.first=grid())

plot(log2norm.prcomp.centred$x[,c(1,4)],
     col=sampleColors,
     main= paste0("PC plot for the experiment no._", parameters$recountID), las=1, panel.first=grid())

par(mfrow=c(1,1))
silence <- dev.off()

#############################################
####### Drax a 3D plot with the 3 first components #####
# install.packages("scatterplot3d", dependencies = TRUE)
library(scatterplot3d)
scatterplot3d.file <- file.path(dir.visualisePCs, paste0("PCA_scatterplot3d_PC3-PC4-PC1_",parameters$recountID,".pdf"))
message("3D scatter plot: for the experiment no.",parameters$recountID," path file :", scatterplot3d.file)
pdf(file = scatterplot3d.file, width=7, height=7)
scatterplot3d(x = log2norm.prcomp.centred$x[,3],
              y = log2norm.prcomp.centred$x[,4],
              z = log2norm.prcomp.centred$x[,1], color=sampleColors,
              xlab="PC3", ylab="PC4", zlab="PC1", las=1, angle=45)
silence <- dev.off()

################################################################
## Another package
# install.packages("rgl", dependencies = TRUE)
# library(rgl)
# plot3d(x = log2norm.prcomp.centred$x[,1],
#               y = log2norm.prcomp.centred$x[,2],
#               z = log2norm.prcomp.centred$x[,3], col = sampleColors,
#               xlab="PC1", ylab="PC2", zlab="PC3", las=1)

## Mustafa, if you want, you can try this (not easy) to get nicer plots.
## https://plot.ly/r/3d-scatter-plots/
# library(plotly)
# p <- plot_ly(as.data.frame(log2norm.prcomp.centred$x[,1:3]))
# chart_link <- plotly_POST(p, filename="scatter3d/basic")
# chart_link

## MUSTAFA: TO CHECK: do you still need the two lines below ? Do you use the variable log2normCounts.components ?
## Seletc the first N components for classification
nb.components <- 4 ## Default number of components
log2normCounts.components <- log2norm.prcomp.centred$x[,1:nb.components]



