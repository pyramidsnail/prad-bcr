library(WGCNA)
library(stringr)
library(stringi)
library(ChAMP)
library(clusterProfiler)
library(org.Hs.eg.db)
library(fgsea)
options(stringsAsFactors = F)

rpkms <- read.table("TCGA-PRAD.htseq_fpkm.tsv.gz", header = T, sep="\t", 
                    row.names = 1)
gencode <- read.table("gencode.v22.annotation.gene.probeMap.gz", header = T)
illuminaMethyl450 <- read.table("illuminaMethyl450_hg38_GDC.gz", header=F)

rpkms <- 2^rpkms-1
code <- substr(colnames(rpkms),14,15) 
rpkms <- rpkms[,as.numeric(code)<10]
filter <- apply(rpkms, 1, mean)
rpkms <- rpkms[which(filter>=1),] 
symbols <- gencode[gencode$id %in% rownames(rpkms), "gene"]
pattern = paste(symbols, collapse="|")
cgs <- illuminaMethyl450[stri_detect_regex(illuminaMethyl450$V2, pattern),]
positions <- gencode[gencode$id %in% rownames(rpkms),]

##### find CpGs which located in the promoter region of genes
cgs.promoter <- cgs[0,]
for (i in seq(1, nrow(positions))){
  cg <- cgs[grepl(positions[i, "gene"],cgs$V2),]
  if (positions[i,"strand"]=="+"){
    start <- as.numeric(positions[i,"chromStart"])
    cg <- cg[as.numeric(cg$V4)<=(start+100) & as.numeric(cg$V4)>=start-2000, ]
  }else{
    start <- as.numeric(positions[i,"chromEnd"])
    cg <- cg[as.numeric(cg$V4)>=(start-100) & as.numeric(cg$V4)<=start+2000, ]
  }
  cgs.promoter <- rbind(cgs.promoter, cg)
}
write.table(cgs.promoter$V1, "data/cg-ids", col.names = F, row.names = F, quote = F)
meths <- read.table("data/illuminaMethyl450-ext", header=F, sep = "\t", row.names = 1)
header <- read.table("data/header", header = F, row.names = 1)
colnames(meths) <- t(header[,3:ncol(header)])
code <- substr(colnames(meths),14,15) 
meths <- meths[,as.numeric(code)<10]
colnames(meths) <- gsub("-",".",colnames(meths))
samples <- intersect(colnames(rpkms), colnames(meths))
rpkms <- rpkms[,colnames(rpkms) %in% samples]
meths <- meths[,colnames(meths) %in% samples]

###### find the negatively correlated pasirs of genes and CpGs
cg.genes <- data.frame(cg=as.character(), gene=as.character(), r=as.character(), p=as.character())
for (i in seq(1, nrow(positions))){
  cg <- cgs[grepl(positions[i, "gene"],cgs$V2),]
  if (positions[i,"strand"]=="+"){
    start <- as.numeric(positions[i,"chromStart"])
    cg <- cg[as.numeric(cg$V4)<=(start+100) & as.numeric(cg$V4)>=start-2000, ]
  }else{
    start <- as.numeric(positions[i,"chromEnd"])
    cg <- cg[as.numeric(cg$V4)>=(start-100) & as.numeric(cg$V4)<=start+2000, ]
  }
  for (j in cg$V1){
    res <- cor.test(t(rpkms[positions[i,"id"],]), t(meths[j, match(colnames(rpkms), colnames(meths))]))
    if (res$p.value<.05 & res$estimate<0){
      cg.genes[nrow(cg.genes)+1,] <-c(j, positions[i,"id"], res$estimate, res$p.value)
      
    }
    
  }
  if (i %% 500==0){
    cat(i)
    cat(" ")
  }
}

##### impute NA values of methylation data
meths.champ <- champ.impute(beta=as.matrix(meths))
meths.champ <- meths.champ$beta
filter <- apply(meths.champ, 1, mean)
filter <- names(filter)[filter>=0.1]
meths.champ <- meths.champ[rownames(meths.champ) %in% filter,]
cg.genes.s <- cg.genes[cg.genes$cg %in% rownames(meths.champ),]
meths.champ <- meths.champ[rownames(meths.champ) %in% cg.genes.s$cg,]

##### aggregate probes on the gene levels 
meths.a <- as.data.frame(meths.champ)
meths.a$id <- cg.genes.s[match(rownames(meths.a), cg.genes.s$cg),"gene"]
meths.a <- aggregate(meths.a, by=list(meths.a$id), median)
rownames(meths.a) <- meths.a$Group.1
meths.a$Group.1 <- NULL
meths.a$id <- NULL
rns <- Reduce(intersect, list(cg.genes.s$gene, rownames(rpkms), rownames(meths.a)))

##### Weighted Gene Coexpression Analysis on the expression levels 
rpkms <- rpkms[rownames(rpkms) %in% rns,]
cor <- WGCNA::cor
datExpr <- log2(t(rpkms+1))
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red")
net.rpkm = blockwiseModules(datExpr, power = 8,maxBlockSize=5000, 
                            TOMType = "unsigned", 
                            reassignThreshold = 0, deepSplit = 4, 
                            numericLabels = T, pamRespectsDendro = FALSE,
                            saveTOMs = TRUE,
                            saveTOMFileBase = "rpkms-tommodel",
                            verbose = 3)
mergedColors = labels2colors(net.rpkm$colors)
plotDendroAndColors(net.rpkm$dendrograms[[1]], mergedColors[net.rpkm$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main="Cluster Dendrogram of Expression")

##### Weighted Gene Coexpression Analysis on the methylation levels 
datExpr <- t(meths.a)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red")
net.meths = blockwiseModules(datExpr, power = 8, maxBlockSize=5000, 
                             TOMType = "unsigned", 
                             reassignThreshold = 0, deepSplit = 4, 
                             numericLabels = T, pamRespectsDendro = FALSE,
                             saveTOMs = TRUE,
                             saveTOMFileBase = "meths-tommodel",
                             verbose = 3)

mergedColors = labels2colors(net.meths$colors)
plotDendroAndColors(net.meths$dendrograms[[1]], mergedColors[net.meths$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main="Cluster Dendrogram of Methylation")

##### Module preservation analysis
multiExpr = list()
multiExpr[[1]] = list(data = log2(t(rpkms+1)))
multiExpr[[2]] = list(data = t(meths.a))
setLabels = c("RNAseq", "methylation")
names(multiExpr) = setLabels
lapply(multiExpr, lapply, dim)
colorList = list(net.rpkm$colors, net.meths$colors)
names(colorList) = setLabels
system.time( {
  mp = modulePreservation(multiExpr, colorList,
                          referenceNetworks = c(1:2),
                          loadPermutedStatistics = FALSE,
                          verbose = 3)
} )

nSets=2
eigengenes = list()
for (set in 1:nSets){
  eigengenes[[set]] = multiSetMEs(multiExpr, universalColors = colorList[[set]], excludeGrey = TRUE);
}
rownames(eigengenes[[1]][[1]]$data) = colnames(rpkms)
rownames(eigengenes[[1]][[2]]$data) = colnames(meths.a)
rownames(eigengenes[[2]][[1]]$data) = colnames(rpkms)
rownames(eigengenes[[2]][[2]]$data) = colnames(meths.a)


### Generate centroids with annotated rows
Centroids <- matrix(rnorm(30, sd = 10), 10)
rownames(Centroids) <- letters[1:nrow(Centroids)]

### Generate data with annotated rows
Data <- cbind(matrix(rep(Centroids[,1], 10), 10),
              matrix(rep(Centroids[,2], 15), 10), matrix(rep(Centroids[,3], 20), 10))
Data <- Data + matrix(rnorm(length(Data), sd = 10), nrow(Data))
rownames(Data) <- letters[1:nrow(Data)]

library(clusterRepro)
cr = list()
set.seed(20)
for (ref in 1:nSets)
{
  cr[[ref]] = list();
  for (test in 1:nSets)
  {
    printFlush(system.time({
      cr[[ref]][[test]] = clusterRepro(Centroids = as.matrix(eigengenes[[ref]][[test]]$data),
                                       New.data = as.matrix(multiExpr[[test]]$data),
                                       Number.of.permutations = 100)}));
    collectGarbage();
  }
}

ref = 1 
test = 2 
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])
print(signif(statsZ[, "Zsummary.pres", drop = FALSE],2));
print(signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))

ref = 1
test = 2
ind = 1
stats= mp$preservation$observed[[ref]][[test]];
labelsX = rownames(stats)
modColors = labelsX;
plotMods = !(modColors %in% c("0", "0.1"));
moduleSizes = stats[plotMods, 1];
textLabels = paste("",modColors, sep="")[plotMods];
colorLabels = labelsX[plotMods];

nModules = sum(plotMods);
nPlots = 6
plotData = list();
# Fill up the plotData
plotData[[1]] = plotData[[2]] = matrix(0, nModules, nPlots);
plotData[[1]][, c(1:4)] = moduleSizes;
plotData[[2]][, 1] = mp$preservation$Z[[ref]][[test]]$Zsummary.pres[plotMods];
plotData[[2]][, 2] = mp$preservation$observed[[ref]][[test]]$medianRank.pres[plotMods];
# Match the modulePreservation ordering of modules to that of clusterRepro
crLabels = sort(unique(colorLabels));
mp2cr = match(colorLabels, crLabels);
# Scatterplots of IGP and p-value vs. module size
plotData[[2]][, 3] = cr[[ref]][[test]]$Actual.IGP
plotData[[2]][, 4] = -log10(cr[[ref]][[test]]$p.value + 1e-4);
# Scatterplot of observed IGP vs. Zsummary and medianRank
plotData[[1]][, c(5,6)] = plotData[[2]][, c(1:2)];
plotData[[2]][, c(5,6)] = plotData[[2]][, 3];
# Plot annotation
xLabs = c(rep("Module size", 4), "Zsummary", "Median rank");
yLabs = c("Zsummary", "Median rank", "Observed IGP", "-log10(IGP perm p)", "Observed IGP", "Observed IGP");
mains = spaste(LETTERS[1:nPlots], ". ", #rep("Ref: Human, Test: Chimp\n", nPlots),
               c(yLabs[1:4], paste(yLabs[5:6], "vs.", xLabs[5:6])),
               c("", "", "", "", "\n", "\n"));
# Scatterplot options
verbose = c(rep(FALSE, 4), rep(TRUE, 2));
ablines = list(c(0, 2, 10), NA, NA, c(-log10(0.05), -log10(0.05/nModules)), NA, NA);
abColors = list(c("black", "blue", "darkgreen"), NA, NA, c("blue", "red"), NA, NA);
logs = c("x", "x", "x", "x", "", "");
invertY = c(FALSE, TRUE, rep(FALSE, 4));
verSP = function(...) { verboseScatterplot(..., abline = TRUE) }

cexLabels = 1.2
sizeGrWindow(6,9);
##### plot Supplemental Figure S2B
pdf(file = "figures/ZStatistics-expr.pdf", w=5, h=5, onefile = FALSE);
for (p in 1){
  x = plotData[[1]][, p];
  y = plotData[[2]][, p]
  miny = min(y, ablines[[p]], na.rm = TRUE);
  maxy = max(y, ablines[[p]], na.rm = TRUE);
  miny = miny - (maxy-miny)*0.1;
  maxy = maxy + (maxy-miny)*0.1;
  (if (verbose[p]) verSP else plot ) (plotData[[1]][, p], plotData[[2]][, p],
                                      xlab = xLabs[p],
                                      ylab = yLabs[p],
                                      cex.main = cexLabels, cex.lab = cexLabels, cex.axis = cexLabels,
                                      bg = "darkgrey",
                                      col = "darkgrey", cex = 2.2,
                                      ylim = c(miny, maxy),
                                      pch = 21,
                                      log = logs[p]);
  labelPoints(plotData[[1]][, p], plotData[[2]][, p], textLabels, cex = cexLabels, offs = 0.00);

  if (!is.na(ablines[[p]][[1]]))
    for (al in 1:length(ablines[[p]]))
      abline(h = ablines[[p]][[al]], col = abColors[[p]][[al]], lty = 2);
}
dev.off()

ref = 2
test = 1
ind = 1
stats= mp$preservation$observed[[ref]][[test]];
labelsX = rownames(stats)
modColors = labelsX;
plotMods = !(modColors %in% c("0", "0.1"));
moduleSizes = stats[plotMods, 1];
textLabels = paste("",modColors, sep="")[plotMods];
colorLabels = labelsX[plotMods];

nModules = sum(plotMods);
nPlots = 6
plotData = list();
# Fill up the plotData
plotData[[1]] = plotData[[2]] = matrix(0, nModules, nPlots);
plotData[[1]][, c(1:4)] = moduleSizes;
plotData[[2]][, 1] = mp$preservation$Z[[ref]][[test]]$Zsummary.pres[plotMods];
plotData[[2]][, 2] = mp$preservation$observed[[ref]][[test]]$medianRank.pres[plotMods];
# Match the modulePreservation ordering of modules to that of clusterRepro
crLabels = sort(unique(colorLabels));
mp2cr = match(colorLabels, crLabels);
# Scatterplots of IGP and p-value vs. module size
plotData[[2]][, 3] = cr[[ref]][[test]]$Actual.IGP
plotData[[2]][, 4] = -log10(cr[[ref]][[test]]$p.value + 1e-4);
# Scatterplot of observed IGP vs. Zsummary and medianRank
plotData[[1]][, c(5,6)] = plotData[[2]][, c(1:2)];
plotData[[2]][, c(5,6)] = plotData[[2]][, 3];
# Plot annotation
xLabs = c(rep("Module size", 4), "Zsummary", "Median rank");
yLabs = c("Zsummary", "Median rank", "Observed IGP", "-log10(IGP perm p)", "Observed IGP", "Observed IGP");
mains = spaste(LETTERS[1:nPlots], ". ", #rep("Ref: Human, Test: Chimp\n", nPlots),
               c(yLabs[1:4], paste(yLabs[5:6], "vs.", xLabs[5:6])),
               c("", "", "", "", "\n", "\n"));
# Scatterplot options
verbose = c(rep(FALSE, 4), rep(TRUE, 2));
ablines = list(c(0, 2, 10), NA, NA, c(-log10(0.05), -log10(0.05/nModules)), NA, NA);
abColors = list(c("black", "blue", "darkgreen"), NA, NA, c("blue", "red"), NA, NA);
logs = c("x", "x", "x", "x", "", "");
invertY = c(FALSE, TRUE, rep(FALSE, 4));
verSP = function(...) { verboseScatterplot(..., abline = TRUE) }

cexLabels = 1.2
sizeGrWindow(6,9);

##### plot Supplemental Figure S2C
pdf(file = "figures/ZStatistics-meths.pdf", w=5, h=5, onefile = FALSE);
for (p in 1){
  x = plotData[[1]][, p];
  y = plotData[[2]][, p]
  miny = min(y, ablines[[p]], na.rm = TRUE);
  maxy = max(y, ablines[[p]], na.rm = TRUE);
  miny = miny - (maxy-miny)*0.1;
  maxy = maxy + (maxy-miny)*0.1;
  (if (verbose[p]) verSP else plot ) (plotData[[1]][, p], plotData[[2]][, p],
                                      xlab = xLabs[p],
                                      ylab = yLabs[p],
                                      cex.main = cexLabels, cex.lab = cexLabels, cex.axis = cexLabels,
                                      bg = "darkgrey",
                                      col = "darkgrey", cex = 2.2,
                                      ylim = c(miny, maxy),
                                      pch = 21,
                                      log = logs[p]);
  labelPoints(plotData[[1]][, p], plotData[[2]][, p], textLabels, cex = cexLabels, offs = 0.00);
  if (!is.na(ablines[[p]][[1]]))
    for (al in 1:length(ablines[[p]]))
      abline(h = ablines[[p]][[al]], col = abColors[[p]][[al]], lty = 2);
}
dev.off()

library(flashClust)
nSets <- 2
dendrograms = list();
inNetwork <- intersect(names(net.meths$colors)[net.meths$colors!=0], names(net.rpkm$colors)[net.rpkm$colors!=0])
for (set in 1:nSets){
  adj = abs(cor(multiExpr[[set]]$data[, inNetwork], use = "p"))^8;
  dtom = TOMdist(adj);
  dendrograms[[set]] = flashClust(as.dist(dtom), method = "a");
}
# Get eigengenes
mes = list()
for (set in 1:nSets)
{
  mes[[set]] = moduleEigengenes(multiExpr[[set]]$data, colorList[[ref]])$eigengenes
}
# Calculate the contingency table and p-values
colorRpkm <- net.rpkm$colors
colorMeths <- net.meths$colors

overlap = overlapTable(colorRpkm[inNetwork], colorMeths[inNetwork]);
# The numMat will encode color. We use -log of the p value.
numMat = -log10(overlap$pTable);
numMat[numMat >50] = 50;
# Prepare for generating a color-coded plot of the overlap table. The text of the table will consist of
# counts and corresponding p-values.
textMat = paste(overlap$countTable, "\n", signif(overlap$pTable, 2));
dim(textMat) = dim(numMat)
# Additional information for the plot. These will be used shortly.
xLabels = paste("cluster", sort(unique(colorMeths)), sep="");
xLabels <- xLabels[2:length(xLabels)]
yLabels = paste("cluster", sort(unique(colorRpkm)), sep="")
yLabels <- yLabels[2:length(yLabels)]

xSymbols = paste(sort(unique(colorMeths)), ": ", table(colorMeths[inNetwork]), sep = "")
xSymbols <- xSymbols[2:length(xSymbols)]
ySymbols = paste(sort(unique(colorRpkm)), ": ", table(colorRpkm[inNetwork]), sep = "")
ySymbols <- ySymbols[2:length(ySymbols)]

##### plot Supplemental Figure S2A
sizeGrWindow(7, 7); fp = FALSE;
pdf("figures/motivationFigure-dendrosAndTable.pdf", w = 10, h = 10)
layout(matrix(c(1,2,5, 3,4,5), 3, 2),
       heights = c(3, 1, 5.5), widths = c(1, 1));
par(mgp = c(3, 1, 0));
plotDendroAndColors(dendrograms[[1]],
                    cbind(colorRpkm[inNetwork], colorMeths[inNetwork]),
                    c("Expression modules", "Methylation modules"),
                    setLayout = FALSE,
                    marAll = c(1, 6, 2.7, 0.2),
                    addGuide = FALSE,
                    main = "A. Expression level dendrogram\nand module colors", cex.main = 1.2,
                    dendroLabels = FALSE, hang = 0.03, cex.colorLabels = 0.7, abHeight = 0.95);

par(mgp = c(3, 1, 0));
plotDendroAndColors(dendrograms[[2]],
                    cbind(colorRpkm[inNetwork], colorMeths[inNetwork]),
                    c("Expression modules", "Methylation modules"),
                    setLayout = FALSE,
                    marAll = c(1, 6, 2.7, 0.2),
                    addGuide = FALSE,
                    main = "B. Methylation level dendrogram\nand module colors", cex.main = 1.2,
                    dendroLabels = FALSE, hang = 0.03, cex.colorLabels = 0.7, abHeight = 0.95);
# Plot the overlap table
fcex = 1.00;
pcex = 1.0
fcexl = 1.00;
pcexl = 1.00;
par(mar = c(6, 7, 2, 1.0));
labeledHeatmap(Matrix = numMat,
               xLabels = xLabels, xSymbols = xSymbols,
               yLabels = yLabels, ySymbols = ySymbols,
               colorLabels = TRUE,
               colors = greenWhiteRed(100)[50:100],
               textMatrix = textMat, cex.text = if (fp) fcex else pcex, setStdMargins = FALSE,
               cex.lab = if (fp) fcexl else pcexl,
               xColorWidth = 0.08,
               main = "C. Expression modules (rows) vs. Methylation modules (columns)", cex.main = 1.2)
dev.off()


##### Correlation analysisi with DNMTs
cnv <- read.table("TCGA-PRAD.gistic.tsv.gz", header=T, row.names = 1, sep="\t")
mutn <- read.table("TCGA-PRAD.mutect2_snv.tsv.gz", header=T,  sep="\t")
rpkms.o <- read.table("TCGA-PRAD.htseq_fpkm.tsv.gz", header = T, sep="\t", 
                      row.names = 1)
rpkms.o <- 2^rpkms.o-1
code <- substr(colnames(rpkms.o),14,15) 
rpkms.o <- rpkms.o[,as.numeric(code)<10]

m5c.genes <- c("DNMT1","DNMT3A","DNMT3B","DNMT3L","MBD1","MBD2","MBD3",'MBD4',"MECP2","NEIL1","NTHL1","SMUG1",
               "TDG","UHRF1","UHRF2","UNG","ZBTB33","ZBTB38","ZBTB4","TET1",'TET2',"TET3")
ids <- gencode[gencode$gene %in% m5c.genes,"id"]
datTraits <- t(rpkms.o[rownames(rpkms.o) %in% ids,])
tempt <- t(cnv[rownames(cnv) %in% ids,])
tempt <- tempt[as.numeric(substr(rownames(tempt),14,15))<10, ]
datTraits <- merge(datTraits, tempt, by=0, all=T)
rownames(datTraits) <- datTraits$Row.names
datTraits$Row.names <- NULL

ensgs <- gsub("(ENSG.*)\\.[x|y]","\\1",colnames(datTraits))
css <- gencode[match(ensgs, gencode$id),"gene"]
css[grep("\\.x$", colnames(datTraits))] <- paste(css[grep("\\.x$", colnames(datTraits))], "exp") 
css[grep("\\.y$", colnames(datTraits))] <- paste(css[grep("\\.y$", colnames(datTraits))], "cnv") 
colnames(datTraits) <- css
datTraits <- datTraits[,c(paste(m5c.genes,"exp"), paste(m5c.genes, "cnv"))]
datTraits <- datTraits[match(colnames(rpkms), rownames(datTraits)),]
datExpr <- log2(t(rpkms)+1)
moduleColors<- net.rpkm$colors
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors,excludeGrey = TRUE)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

textSymbol <- textMatrix
textSymbol[moduleTraitPvalue<0.01] <- "*"
textSymbol[moduleTraitPvalue>0.01 & moduleTraitPvalue<.05] <- "+"
textSymbol[moduleTraitPvalue>.05] <- ""

# Display the correlation values within a heatmap plot, Figure3
pdf("figures/expr-genes-trait.pdf", height = 5, width=12)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = gsub("ME","cluster",names(MEs)),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textSymbol,
               setStdMargins = FALSE,
               cex.text = 1.1,
               zlim = c(-1,1),
               main = "")
dev.off()

##### for the methylation clusters
datTraits <- datTraits[match(colnames(meths.a), rownames(datTraits)),]
datExpr <- t(meths.a)
moduleColors<- net.meths$colors
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
MEs0 = moduleEigengenes(datExpr, moduleColors,excludeGrey = TRUE)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

textSymbol <- textMatrix
textSymbol[moduleTraitPvalue<0.01] <- "*"
textSymbol[moduleTraitPvalue>0.01 & moduleTraitPvalue<.05] <- "+"
textSymbol[moduleTraitPvalue>.05] <- ""
# Display the correlation values within a heatmap plot, Supplemental Figure S1
pdf("figures/meths-genes-trait.pdf", height = 4, width=12)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = gsub("ME","cluster",names(MEs)),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textSymbol,
               setStdMargins = FALSE,
               cex.text = 1.1,
               zlim = c(-1,1),
               main = "")
dev.off()

##### functional enrichment analysis for each coexpression and comethylation module 
##### generating Figure 2B, 2D, Figure 4A, 4B, also supplemental tables
library(clusterProfiler)
library(org.Hs.eg.db)
hallmark <- gmtPathways("h.all.v7.0.symbols.gmt")

datExpr <- log2(t(rpkms)+1)
moduleColors<- net.rpkm$colors
MEs0 = moduleEigengenes(datExpr, moduleColors,excludeGrey = TRUE)$eigengenes
MEs = orderMEs(MEs0)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
for ( i in seq(1, 19)){
  genes  <- names(net.rpkm$colors)[net.rpkm$colors==i]
  genes <- gencode[gencode$id %in% genes, "gene"]
  go <- enrichGO(genes, OrgDb = "org.Hs.eg.db", ont="all", keyType = "SYMBOL")
  
  genes <- names(net.rpkm$colors)[net.rpkm$colors==i]
  colname <- paste("ME", i, sep="")
  ranks <- geneModuleMembership[genes, colname]
  names(ranks) <- gencode[match(genes, gencode$id), "gene"]
  
  hm <- fgsea(hallmark,ranks,nperm=1000,minSize=5)
  assign(paste("expr",i,"hm", sep="."), hm)
  assign(paste("expr",i,"go", sep="."), go)
  
}


datExpr <- t(meths.a)
moduleColors<- net.meths$colors
MEs0 = moduleEigengenes(datExpr, moduleColors,excludeGrey = TRUE)$eigengenes
MEs = orderMEs(MEs0)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

for ( i in seq(1, 13)){
  genes  <- names(net.meths$colors)[net.meths$colors==i]
  genes <- gencode[gencode$id %in% genes, "gene"]
  go <- enrichGO(genes, OrgDb = "org.Hs.eg.db", ont="all", keyType = "SYMBOL")
  
  genes <- names(net.meths$colors)[net.meths$colors==i]
  colname <- paste("ME", i, sep="")
  ranks <- geneModuleMembership[genes, colname]
  names(ranks) <- gencode[match(genes, gencode$id), "gene"]
  
  hm <- fgsea(hallmark,ranks,nperm=1000,minSize=5)
  assign(paste("meths",i,"hm", sep="."), hm)
  assign(paste("meths",i,"go", sep="."), go)
}

##### calculate module-trait correlation for the selected traits 
survival <- read.table("~/projects/golden/prad/data/TCGA-PRAD.survival.tsv.gz", header = T, sep="\t")
phenotype <- read.table("~/projects/golden/prad/data/TCGA-PRAD.GDC_phenotype.tsv.gz", header = T, sep="\t", quote = "\"")

datTraits <- phenotype[,c("submitter_id.samples","pathologic_T","pathologic_N", "clinical_T","biochemical_recurrence",
                          "days_to_first_biochemical_recurrence", "days_to_second_biochemical_recurrence","days_to_third_biochemical_recurrence",
                          "gleason_score","psa_value")]
datTraits$pathologic_T <- gsub(".*([0-9]).*", "\\1", datTraits$pathologic_T)
datTraits$pathologic_N <- gsub(".*([0-9]).*", "\\1", datTraits$pathologic_N)
datTraits$clinical_T <- gsub(".*([0-9]).*", "\\1", datTraits$clinical_T)
datTraits$biochemical_recurrence[datTraits$biochemical_recurrence=="YES"] <-1
datTraits$biochemical_recurrence[datTraits$biochemical_recurrence=="NO"] <-0
datTraits$submitter_id.samples <- gsub("-",".", datTraits$submitter_id.samples)
rownames(datTraits) <- datTraits$submitter_id.samples
datTraits$submitter_id.samples <- NULL
datTraits <- datTraits[match(colnames(rpkms), rownames(datTraits)),]
datExpr <- log2(t(rpkms)+1)
moduleColors<- net.rpkm$colors
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
MEs0 = moduleEigengenes(datExpr, moduleColors,excludeGrey = TRUE)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

#####  generating Figure 2A
pdf("figures/expr-trait.pdf", height = 7, width=7)
par(mar=c(8,5,2,2))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = gsub("_biochemical_recurrence","\nbiochemical_recurrence",names(datTraits)),
               yLabels = gsub("ME","cluster",names(MEs)),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = .6,
               zlim = c(-1,1),
               main = "")
dev.off()


datTraits <- datTraits[match(colnames(meths.a), rownames(datTraits)),]
datExpr <- t(meths.a)
moduleColors<- net.meths$colors
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
MEs0 = moduleEigengenes(datExpr, moduleColors,excludeGrey = TRUE)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

#####  generating Figure 2C
pdf("figures/meths-trait.pdf", height = 7, width=7)
par(mar=c(8,5,2,2))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = gsub("_biochemical_recurrence","\nbiochemical_recurrence",names(datTraits)),
               yLabels = gsub("ME","cluster",names(MEs)),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = .6,
               zlim = c(-1,1),
               main = "")
dev.off()

##### for immune cell infiltration analysis
colnames(expr.g)[1] <- "Gene"
write.table(expr.g, "PRAD-mixture.txt", col.names=T,row.names=F, sep="\t", quote=F)
source("CIBERSORT.R")
results <- CIBERSORT("LM22.txt","PAAD-mixture.txt", perm=100)
cs <- results
datTraits <- cs[,1:22]

datTraits <- datTraits[match(colnames(rpkms), rownames(datTraits)),]
MEs0 = moduleEigengenes(datExpr.rpkm, moduleColors.rpkm,excludeGrey = TRUE)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

textSymbol <- textMatrix
textSymbol[moduleTraitPvalue<0.01] <- "*"
textSymbol[moduleTraitPvalue>0.01 & moduleTraitPvalue<.05] <- "+"
textSymbol[moduleTraitPvalue>.05] <- ""

#####  generating Supplemental Figure S3A
pdf("figures/expr-CIBERSORT.pdf", height = 6, width=8)
par(mar=c(10,5,1,2))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits),
               yLabels = gsub("ME","cluster",names(MEs)),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textSymbol,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = "")
dev.off()

datTraits <- datTraits[match(colnames(meths.a), rownames(datTraits)),]
MEs0 = moduleEigengenes(datExpr.meth, moduleColors.meth,excludeGrey = TRUE)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
textSymbol <- textMatrix
textSymbol[moduleTraitPvalue<0.01] <- "*"
textSymbol[moduleTraitPvalue>0.01 & moduleTraitPvalue<.05] <- "+"
textSymbol[moduleTraitPvalue>.05] <- ""

#####  #####  generating Supplemental Figure S3B
pdf("figures/meths-CIBERSORT.pdf", height = 5, width=8)
par(mar=c(10,5,1,2))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits),
               yLabels = gsub("ME","cluster",names(MEs)),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textSymbol,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = "")
dev.off()

##### TILs analysis with 5mC regulator expressions
ids <- reference[reference$gene %in% m5c.genes, "id"]
rpkms.dnmt <- t(rpkms.o[ids,])
dnmt.til.cor <- data.frame(til=as.character(), dnmt=as.character(), coef=as.numeric(), p=as.numeric())
for (i in seq(1,22)){
  for (j in colnames(rpkms.dnmt)){
    test <- cor.test(cs[,i], rpkms.dnmt[match(rownames(cs), rownames(rpkms.dnmt)),j])
    symbol <- reference[reference$id==j, "gene"]
    dnmt.til.cor[nrow(dnmt.til.cor)+1,] <- c(colnames(cs)[i], symbol, test$estimate, test$p.value)
  }
}
dnmt.til.cor$symbol <- ""
dnmt.til.cor$symbol[as.numeric(dnmt.til.cor$p)<0.1] <- "+"
dnmt.til.cor$symbol[as.numeric(dnmt.til.cor$p)<0.05] <- "*"
dnmt.til.cor$coef <- as.numeric(dnmt.til.cor$coef)
dnmt.til.cor$dnmt <- factor(dnmt.til.cor$dnmt, levels = rev(m5c.genes))
p <- ggplot(dnmt.til.cor, aes(til, dnmt, fill=coef))+geom_tile()+
  scale_fill_gradient2(low="steelblue",mid="white", high="darkred",midpoint = 0, name="correlation\ncoefficient")+
  theme(axis.text.x = element_text(angle = 90, size=18, hjust = 1, vjust = .5), axis.text.y = element_text(size=18), legend.title = element_text(size=18), legend.text = element_text(size=15))+
  geom_text(aes(label=symbol),size=7)+xlab("")+ylab("")

#####  generating Figure 6A
ggsave("figures/dnmt-tils.pdf", height=10, width=10)
dnmt.til.cor[abs(as.numeric(dnmt.til.cor$coef))>.3 & as.numeric(dnmt.til.cor$p)<.05,]

##### CPBs with 5mC regulator expressions
cpbs <- c("PDCD1","CD274","LAG3", "HAVCR2", "CD47", "IDO1","CTLA4","CD80","VTCN1","ADORA2A","ARHGEF5","BTLA","CD160","CD244","CD27","CD276",
          "CEACAM1","GEM","ICOS","TNFSF4","VISTA")
cpbs.ids <- reference[reference$gene %in% cpbs, "id"]
ids <- reference[reference$gene %in% m5c.genes, "id"]
rpkms.dnmt <- t(rpkms.o)
dnmt.cpb.cor <- data.frame(cpb=as.character(), dnmt=as.character(), coef=as.numeric(), p=as.numeric())
for (i in cpbs.ids){
  for (j in ids){
    test <- cor.test(rpkms.dnmt[,i], rpkms.dnmt[,j])
    symbol1 <- reference[reference$id==i, "gene"]
    symbol2 <- reference[reference$id==j, "gene"]
    dnmt.cpb.cor[nrow(dnmt.cpb.cor)+1,] <- c(symbol1, symbol2, test$estimate, test$p.value)
  }
}
dnmt.cpb.cor$symbol <- ""
dnmt.cpb.cor$symbol[as.numeric(dnmt.cpb.cor$p)<0.1] <- "+"
dnmt.cpb.cor$symbol[as.numeric(dnmt.cpb.cor$p)<0.05] <- "*"
dnmt.cpb.cor$coef <- as.numeric(dnmt.cpb.cor$coef)
dnmt.cpb.cor$dnmt <- factor(dnmt.cpb.cor$dnmt, levels = rev(m5c.genes))
p <- ggplot(dnmt.cpb.cor, aes(cpb, dnmt, fill=coef))+geom_tile()+
  scale_fill_gradient2(low="steelblue",mid="white", high="darkred",midpoint = 0, name="correlation\ncoefficient")+
  theme(axis.text.x = element_text(angle = 90, size=18, hjust = 1, vjust = .5), axis.text.y = element_text(size=18), legend.title = element_text(size=18), legend.text = element_text(size=15))+
  geom_text(aes(label=symbol),size=7)+xlab("")+ylab("")

#####  generating Figure 6B
ggsave("figures/dnmt-cpbs.pdf", height=10, width=10)
dnmt.cpb.cor[abs(as.numeric(dnmt.cpb.cor$coef))>.3 & as.numeric(dnmt.cpb.cor$p)<.05,]


##### dnmts with pathologic phenotypes 
datTraits <- phenotype[,c("submitter_id.samples","pathologic_T","pathologic_N", "clinical_T","biochemical_recurrence",
                          "days_to_first_biochemical_recurrence", "days_to_second_biochemical_recurrence","days_to_third_biochemical_recurrence",
                          "gleason_score","psa_value")]
datTraits$pathologic_T <- gsub(".*([0-9]).*", "\\1", datTraits$pathologic_T)
datTraits$pathologic_N <- gsub(".*([0-9]).*", "\\1", datTraits$pathologic_N)
datTraits$clinical_T <- gsub(".*([0-9]).*", "\\1", datTraits$clinical_T)
datTraits$biochemical_recurrence[datTraits$biochemical_recurrence=="YES"] <-1
datTraits$biochemical_recurrence[datTraits$biochemical_recurrence=="NO"] <-0
datTraits$submitter_id.samples <- gsub("-",".", datTraits$submitter_id.samples)
rownames(datTraits) <- datTraits$submitter_id.samples
datTraits$submitter_id.samples <- NULL
ids <- reference[reference$gene %in% m5c.genes, "id"]
rpkms.dnmt <- t(rpkms.o[ids,])

dnmt.pheno.cor <- data.frame(phenotype=as.character(), dnmt=as.character(), coef=as.numeric(), p=as.numeric())
for (i in colnames(datTraits)){
  for (j in colnames(rpkms.dnmt)){
    test <- cor.test(as.numeric(as.character(datTraits[,i])), rpkms.dnmt[match(rownames(datTraits), 
                                                                               rownames(rpkms.dnmt)),j])
    symbol <- reference[reference$id==j, "gene"]
    dnmt.pheno.cor[nrow(dnmt.pheno.cor)+1,] <- c(i, symbol, test$estimate, test$p.value)
  }
}
dnmt.pheno.cor$symbol <- ""
dnmt.pheno.cor$symbol[as.numeric(dnmt.pheno.cor$p)<0.1] <- "+"
dnmt.pheno.cor$symbol[as.numeric(dnmt.pheno.cor$p)<0.05] <- "*"
dnmt.pheno.cor$coef <- as.numeric(dnmt.pheno.cor$coef)
dnmt.pheno.cor$coef[is.na(dnmt.pheno.cor$coef)] <-0
dnmt.pheno.cor$phenotype <- gsub("_biochemical_recurrence","\nbiochemical_recurrence",dnmt.pheno.cor$phenotype)
dnmt.pheno.cor$phenotype <- factor(dnmt.pheno.cor$phenotype, levels = gsub("_biochemical_recurrence","\nbiochemical_recurrence",colnames(datTraits)))
dnmt.pheno.cor$dnmt <- factor(dnmt.pheno.cor$dnmt, levels = rev(m5c.genes))
p <- ggplot(dnmt.pheno.cor, aes(phenotype, dnmt, fill=coef))+geom_tile()+theme_bw()+
  scale_fill_gradient2(low="steelblue",mid="white", high="darkred",midpoint = 0, name="correlation\ncoefficient")+
  theme(axis.text.x = element_text(angle = 90, size=24, hjust = 1, vjust = .5), axis.text.y = element_text(size=18), legend.title = element_text(size=18), legend.text = element_text(size=15))+
  geom_text(aes(label=symbol),size=7)+xlab("")+ylab("")

#####  generating Figure 6C
ggsave("figures/dnmt-phenotype.pdf", height=10, width=10)
dnmt.pheno.cor[abs(as.numeric(dnmt.pheno.cor$coef))>.3 & as.numeric(dnmt.pheno.cor$p)<.05,]

##### analysis 5mC regulators with BCR
traits <- phenotype
subs <- (traits$biochemical_recurrence=="NO"| traits$biochemical_recurrence=="") & as.numeric(traits$days_to_last_follow_up.diagnoses)<3*365
traits <- traits[-which(subs==TRUE),]
subs <- traits$biochemical_recurrence=="" & is.na(traits$days_to_first_biochemical_recurrence)
traits <- traits[-which(subs==TRUE),]
subs <- traits$biochemical_recurrence=="YES" & is.na(traits$days_to_first_biochemical_recurrence)
traits <- traits[-which(subs==TRUE),]
traits$status <- 0
traits$status[!is.na(traits$days_to_first_biochemical_recurrence)] <- 1
traits$time <- traits$days_to_first_biochemical_recurrence
traits$time[is.na(traits$time)] <- traits[is.na(traits$time),"days_to_last_follow_up.diagnoses"]
traits <- traits[complete.cases(traits[,c("time", "status")]),]

library(survminer)
library(survival)
library(pROC)
library("magrittr")
library("dplyr")
samples <- intersect(gsub("-",".",traits$submitter_id.samples), colnames(rpkms.o))
traits$submitter_id.samples <- gsub("-",".",traits$submitter_id.samples)
traits <- traits[traits$submitter_id.samples %in% samples,]

##### univariate cox regression analysis
cox.dnmt <- data.frame(gene = as.character(), coef = as.character(), p=as.character())
for (id in ids){
  tmpt <- cbind(traits[,c("time","status")], t(rpkms.o[id, match(traits$submitter_id.samples, colnames(rpkms.o))]))
  colnames(tmpt)[ncol(tmpt)] <- "dnmt"
  res.cox <- summary(coxph(Surv(time, status) ~ dnmt, data = tmpt))
  symbol <- reference[reference$id==id,"gene"]
  cox.dnmt[nrow(cox.dnmt)+1,] <- c(symbol, res.cox$coefficients[2],res.cox$coefficients[5])
  
}

genes <- cox.dnmt[as.numeric(cox.dnmt$p)<0.05,"gene"]
ids <- reference[reference$gene %in% genes, "id"]
tmpt <- cbind(traits[,c("time","status")], t(rpkms.o[ids, match(traits$submitter_id.samples, colnames(rpkms.o))]))
colnames(tmpt)[3:ncol(tmpt)] <- reference[match(colnames(tmpt)[3:ncol(tmpt)], reference$id), "gene"]

cox.dnmt.utils <- cox.res.df[0,]
endpoint <- "time"
endpoint.code <- "status"

for (id in ids){
  tmpt <- cbind(traits[,c("time","status")], t(rpkms.o[id, match(traits$submitter_id.samples, colnames(rpkms.o))]))
  symbol <- reference[reference$id==id,"gene"]
  colnames(tmpt)[ncol(tmpt)] <- symbol
  feature <- symbol
  cox.res.df <- get_cox_res(tmpt, endpoint, endpoint.code, feature )
  cox.dnmt.utils[nrow(cox.dnmt.utils)+1,] <- cox.res.df
  
}

plot_cox_res(cox.dnmt.utils[cox.dnmt.utils$p.value<0.05,],x.lab = "Hazard Ratio", 
             y.lab = "")+theme_bw()+theme(axis.text = element_text(size=11, face="bold", color="black"))

##### Figure 6A
ggsave("figures/cox-univariate.pdf", height = 3.5, width=4)
write.table(cox.dnmt.utils, "tables/cox-univariate.csv", col.names = T, row.names = F, quote=F, sep=",")

##### build BCR prediction model 
library(caret)
library(glmnet)
sample.i <-createDataPartition(tmpt$status, p = .7, list = FALSE, times = 1)
train.set <- tmpt[sample.i,]
test.set <- tmpt[-sample.i,]
cv_output <- cv.glmnet(as.matrix(train.set[,3:(ncol(train.set))]),Surv(train.set$time, train.set$status), 
                       family="cox")
plot(cv_output)

fit = glmnet(as.matrix(train.set[,3:(ncol(train.set))]),Surv(train.set$time, train.set$status), 
             family="cox")
plot(fit)

##### get model coef
coef(cv_output, s = "lambda.min")
rs.df <- train.set
##### coef get from coef(cv_output, s = "lambda.min")
rs.score <- 0.48*rs.df$DNMT3B+0.04*rs.df$DNMT1

rs.df$score <- rs.score
rs.df$sample <- traits$submitter_id.samples[sample.i]
rs.df <- cbind(rs.df, datTraits[match(rs.df$sample, rownames(datTraits)),])
rs.df$age <- phenotype[match(rs.df$sample, gsub("-",".", phenotype$submitter_id.samples)), "age_at_diagnosis.diagnoses"]
rs.df$age <- round(as.numeric(rs.df$age)/365,0)


dw <- rs.df[,c( "status", "age","score","gleason_score","psa_value")]
dw <- melt(dw, id.vars = c("status"))
dw$status[dw$status==0] <- "NO"
dw$status[dw$status==1] <- "YES"

dummy <- data.frame(status = "YES", value = 10,
                    variable = "gleason_score", stringsAsFactors=FALSE)
p <- ggplot(dw, aes(factor(status), value))+geom_boxplot()+geom_jitter(alpha=.4)+stat_compare_means(label.y.npc="top")+
  geom_blank(data=dummy) + theme_bw()+
  facet_wrap(.~variable, scales = "free",nrow=1)+ylab("")+xlab("Biochemical Recurrence")
##### figure 7D
ggsave("figures/train-boxplot.pdf", height=2.5, width=9)

rocobj <- roc(rs.df$status, rs.df$score)
g <- ggroc(rocobj,size = 1)+theme_bw()+annotate("text",x=0.5,y=0.1, label="AUC=0.70",size=6)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed")
##### figure 6C
ggsave("figures/train-roc.pdf", height=3, width=3)

roc(rs.df$status, rs.df$gleason_score)
roc(rs.df$status, rs.df$psa_value)
res.cox <- coxph(Surv(time, status) ~ score+gleason_score+psa_value+age, data = rs.df)
ggforest(res.cox) 
##### figure 7B
ggsave("figures/train-cox-forest.pdf", height=3, width=6)

res.cox <- coxph(Surv(time, status) ~ gleason_score+psa_value+age, data = rs.df)
ggforest(res.cox) 
###### figure 7C
ggsave("figures/train-cox-forest-without-score.pdf", height=2, width=6)


cox.train.utils <- cox.res.df[0,]
endpoint <- "time"
endpoint.code <- "status"

for (id in c("score","gleason_score","psa_value","age")){
  
  cox.res.df <- get_cox_res(rs.df, endpoint, endpoint.code, id )
  cox.train.utils[nrow(cox.train.utils)+1,] <- cox.res.df
  
}
write.table(cox.train.utils, "tables/cox-train-utils.csv", col.names = T, row.names = F, quote=F, sep=",")

plot_cox_res(cox.train.utils,x.lab = "Hazard Ratio", 
             y.lab = "")+theme_bw()+theme(axis.text = element_text(size=11, face="bold", color="black"))
##### figure 7A
ggsave("figures/train-cox-univariate.pdf", height = 2, width=4)

mid <- quantile(rs.df$score, seq(0,1,length.out = 3))[2]
rs.df$category <- "Low"
rs.df$category[rs.df$score>mid] <- "High"
fit <- survfit(Surv(time, status) ~ category, data=rs.df)
p.value <- surv_pvalue(fit)$pval
p <- ggsurvplot(fit, conf.int = T, pval=T,
                ggtheme = theme_minimal()+theme(plot.title = element_text(hjust = 0.5, face = "bold")),                   
                legend="right", legend.title="",
                legend.labs = c("High","Low") ,
                palette =pal_npg("nrc", alpha = 0.7)(2),risk.table = TRUE)
##### figure 6B
pdf("figures/train-pfs.pdf", height = 4.5, width=4, onefile = F)
print(print(p))
dev.off()

##### testing on the validaton datasets
rs.df <- test.set
rs.score <- 0.48*rs.df$DNMT3B+0.04*rs.df$DNMT1

rs.df$score <- rs.score
rs.df$sample <- traits$submitter_id.samples[-sample.i]
rs.df <- cbind(rs.df, datTraits[match(rs.df$sample, rownames(datTraits)),])
rs.df$age <- phenotype[match(rs.df$sample, gsub("-",".", phenotype$submitter_id.samples)), "age_at_diagnosis.diagnoses"]
rs.df$age <- round(as.numeric(rs.df$age)/365,0)

dw <- rs.df[,c( "status", "age","score","gleason_score","psa_value")]
dw <- melt(dw, id.vars = c("status"))
dw$status[dw$status==0] <- "NO"
dw$status[dw$status==1] <- "YES"

dummy <- data.frame(status = "YES", value = 10,
                    variable = "gleason_score", stringsAsFactors=FALSE)
p <- ggplot(dw, aes(factor(status), value))+geom_boxplot()+geom_jitter(alpha=.4)+stat_compare_means(label.y.npc="top")+
  geom_blank(data=dummy) + theme_bw()+
  facet_wrap(.~variable, scales = "free",nrow=1)+ylab("")+xlab("Biochemical Recurrence")
#####  figure 8D
ggsave("figures/test-boxplot.pdf", height=2.5, width=9)


rocobj <- roc(rs.df$status, rs.df$score)
roc(rs.df$status, rs.df$score)
g <- ggroc(rocobj,size = 1)+theme_bw()+annotate("text",x=0.5,y=0.1, label="AUC=0.88",size=6)+
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed")
#####  figure 6E
ggsave("figures/test-roc.pdf", height=3, width=3)


roc(rs.df$status, rs.df$gleason_score)
roc(rs.df$status, rs.df$psa_value)
res.cox <- coxph(Surv(time, status) ~ score+gleason_score+psa_value+age, data = rs.df)
ggforest(res.cox) 
#####  figure 8B
ggsave("figures/test-cox-forest.pdf", height=3, width=6)

res.cox <- coxph(Surv(time, status) ~ gleason_score+psa_value+age, data = rs.df)
ggforest(res.cox) 
#####  figure 8C
ggsave("figures/test-cox-forest-withoutscore.pdf", height=2, width=6)

cox.test.utils <- cox.res.df[0,]
endpoint <- "time"
endpoint.code <- "status"

for (id in c("score","gleason_score","psa_value","age")){
  
  cox.res.df <- get_cox_res(rs.df, endpoint, endpoint.code, id )
  cox.test.utils[nrow(cox.test.utils)+1,] <- cox.res.df
  
}
write.table(cox.test.utils, "tables/cox-test-utils.csv", col.names = T, row.names = F, quote=F, sep=",")

plot_cox_res(cox.test.utils,x.lab = "Hazard Ratio", 
             y.lab = "")+theme_bw()+theme(axis.text = element_text(size=11, face="bold", color="black"))
#####  figure 8A
ggsave("figures/test-cox-univariate.pdf", height = 2, width=4)


mid <- quantile(rs.df$score, seq(0,1,length.out = 3))[2]
rs.df$category <- "Low"
rs.df$category[rs.df$score>mid] <- "High"
fit <- survfit(Surv(time, status) ~ category, data=rs.df)
p.value <- surv_pvalue(fit)$pval
p <- ggsurvplot(fit, conf.int = T, pval=T,
                ggtheme = theme_minimal()+theme(plot.title = element_text(hjust = 0.5, face = "bold")),                   
                legend="right", legend.title="",
                legend.labs = c("High","Low") ,
                palette =pal_npg("nrc", alpha = 0.7)(2),risk.table = TRUE)
#####  figure 6D
pdf("figures/test-pfs.pdf", height = 4.5, width=4, onefile = F)
print(print(p))
dev.off()