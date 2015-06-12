#####Load Data //Required
- library("limma","lumi","AnnotationDbi","lumiMouseAll.db","rgl","annotate","R2HTML")
- setwd("/home/alex/Documents/**Illumina**")
- Make /data/ directory with .CEL files
- Make Phenod.txt file with tab-delimited text in Excel with columns Sample,Species,Type,Treatment,Time,Replicate and one row for each dataset

#####Normalize Data //Required
- x.lumi <- lumiR('**GSE43221_non-normalized_data.txt**', lib='lumiMouseAll')
- phenod<-read.table("Phenod.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1, sep="\t")
- lumi.T <- lumiT(x.lumi)
- lumi.N <- lumiN(lumi.T, method='rsn')
- lumi.N.Q <- lumiQ(lumi.N)

######Quality Control Normalized Data
- boxplot(x.lumi)
- boxplot(lumi.N.Q)
- hist(x.lumi)
- hist(lumi.N.Q)

######PCA
- pca<-prcomp(t(x.lumi),scale=T)
- timeseries.colors <- c(rep("Green", 3),
rep("DimGrey", 3),
rep("IndianRed",3),
rep("Red", 3), #WT samples
rep("BurlyWood", 3),
rep("grey", 3),
rep("pink", 3),
rep("SkyBlue", 3)#KO samples
)
- plot3d(pca$x[, 1:3], col=timeseries.colors, xlab="PC1", ylab = "PC2", zlab = "PC3", type = "s")
- rgl.postscript("PCA.pdf", fmt="pdf", drawText=TRUE)

#####Filter Data //Required
- lumi.N.Q1 <- lumiExpresso(x.lumi, normalize.param=list(method='rsn'))
- dataMatrix <- exprs(lumi.N.Q)
- write.table(dataMatrix, file="Data-Matrix-NQ.txt", quote= FALSE, sep = "\t")
- presentCount <- detectionCall(x.lumi)
- selDataMatrix <- dataMatrix[presentCount > **19**,]
- write.table(selDataMatrix,file="Data-Matrix-NQ-filtered.txt", quote=FALSE, sep = "\t")
- probeList <- rownames(selDataMatrix)

#####Select Two Sets (Six Rows) of Data //Required
- mData <- selDataMatrix[, **c(1,2,3,10,11,12)**]
- mDataPheno<-phenod[**c(1,2,3,10,11,12)**, 4]
- design <- model.matrix(~ factor(mDataPheno, levels=c("**Ctrl**", "**24hr**")))
- colnames(design) <- c("**Ctrl**", "**24hrvsCtrl**")

#####Generate Linear Model //Required
- fit <- lmFit(mData, design)
- fit <- eBayes(fit)
- tab <- topTable(fit, coef = 2, adjust = "fdr", n = dim(mData)[1])
- ID <- rownames(tab)
- tab1 <- cbind(ID, tab)
- write.table(tab1, file="Allprobes-pvalues.xls", row.names=FALSE, quote=FALSE, sep="\t")

#####Select Probesets //Required
- tab.sig <- tab[tab$adj.P.Val < **0.05** & tab$logFC > abs(log2(**2**)),]
- ID <- rownames(tab.sig)
- tab.sig1 <- cbind(ID, tab.sig)
- write.table(tab.sig1, file="DEG-fdr0.001.xls", row.names=FALSE, quote=FALSE, sep="\t")

######Generate Heatmap
- selected  <- rownames(tab.sig)
- mDataSig<-mData[rownames(mData) %in% selected,]
- pdf("Clustering.pdf")
- heatmap(mDataSig, labRow=c(""), col=topo.colors(16), cexCol=0.6)
- dev.off()

#####Annotate Results Output //Required
- nuID <- rownames(tab.sig)
- head(nuID)
- Symbol <- getSYMBOL(nuID, 'lumiMouseAll')
- Name <- as.character(lookUp(nuID, 'lumiMouseAll', "GENENAME"))
- Ensembl <- as.character(lookUp(nuID, 'lumiMouseAll', "ENSEMBL"))
- Ensembl <- ifelse(Ensembl=="NA", NA, paste("<a **href** "useast.ensembl.org/Mus_musculus/Gene/Summary?g=", Ensembl, "'>", Ensembl, "</a>", sep=""))
- tmp <- data.frame(ID=nuID, Symbol=Symbol, Name=Name, logFC=tab.sig$logFC, pValue=tab.sig$P.Value, FDR=tab.sig$adj.P.Val, Ensembl=Ensembl)
- row.names(tmp)<-NULL
- HTML(tmp, "out.html", append=FALSE)
- write.table(tmp, file="target.txt", quote=FALSE, row.names=FALSE, sep="\t")

######Clean Up Workspace
- save.image()
- sessionInfo()
