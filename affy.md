#####Load Data //Required
- library("GEOquery","simpleaffy","RColorBrewer","affyPLM","rgl","limma","annotate","hgu133plus2.db")
- setwd("/home/alex/Documents/**HUVEC**")
- Make /data/ directory with .CEL files
- Make phenodata.txt file with tab-delimited text in Excel with columns Name,Filename,Target and one row for each .CEL file

#####Normalize Data //Required
- cellfiles<-read.affy(covdesc="phenodata.txt",path="data")
- cellfiles.gcrma<-gcrma(cellfiles)

#####Quality Control Normalized Data
- boxplot(cellfiles,col=rainbow(12))
- boxplot(cellfiles.gcrma,col=rainbow(12))
- hist(cellfiles,col=rainbow(12))
- hist(cellfiles.gcrma,col=rainbow(12))

#####Clusters
- eset<-exprs(celfiles.gcrma)
- distance<-dist(t(eset),method="maximum")
- clusters<-hclust(distance)
- plot(clusters)

#####Filter Data //Required
- cellfiles.filtered<-nsFilter(cellfiles.gcrma,require.entrez=F,remove.dupEntrez=F)
- filterEset<-exprs(cellfiles.filtered$eset)

#####Generate Linear Models //Required
- samples<-as.factor(cellfiles.gcrma$Target)
- design<-model.matrix(~0+samples)
- levels.design<-levels(samples)
- colnames(design)<-levels.design
- fit<-lmFit(filterEset,design)
- contrast.matrix<-makeContrasts(contrast=**huvec**-**choroid**,levels=design)
- fits<-contrast.fit(fit,contrast.matrix)
- ebFits<-eBayes(fits)
- probeset.list.general<-topTable(ebFits,coef=1,number=nrow(filterEset),lfc=0)

#####Select Probesets //Required
- probeset.list.selected<-probeset.list.general[(probeset.list$adj.P.Val<=**0.05**)&(abs(probeset.list.general$logFC)>=**2**), ]
- cellSamples<-subset(cellfiles.gcrma$FileName,(cellfiles.gcrma$Target=="**huvec**")|(cellfiles.gcrma$Target=="**choroid**"))
- cellData<-filterEset[rownames(filterEset) %in% rownames(probeset.list.selected),colnames(filterEset) %in% cellSamples]

#####Generate Heatmap
- pdf(file="Heatmap.pdf",width=8,height=6)
- heatmap(cellData,labRow=c(""),col=topo.colors(16),cexCol=0.6)
- graphics.off()

#####Annotate Results Output //Required
- gene.symbols<-getSYMBOL(rownames(probeset.list.selected), "hgu133plus2")
- results<-cbind(probeset.list.selected,gene.symbols)
- write.table(results,"results.txt",sep="\t",quote=F)

#####Clean Up Workspace
- save.image()
- sessionInfo()
