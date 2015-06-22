geneProfilingHumans=function(input.expressionfile,sep,input.phenofile,norm,filter.bool,filter.int,phenoline.start,phenoline.finish,phenoname.ctrl,phenoname.var,pval,fcval,heatmap.out,genes.out){
  require("limma")
  require("lumi")
  require("AnnotationDbi")
  require("lumiMouseAll.db")
  require("rgl")
  require("ggplot2")
  require("annotate")
  require("R2HTML")
  #message ("Example Parameters: geneProfilingIllumina("S24-GSdata-noNorm.txt",",","Phenod.txt","rsn",F,19,1,6,"CEM_PLAIN","CEM_SUV",0.20,1.3,"Heatmap-test.pdf","Genes-test.txt")")
  phenod<-read.table(input.phenofile, header=T, stringsAsFactors=FALSE, row.names=1, sep=sep)
  message ("Read Pheno Data")
  annotation=read.table("annotation.txt",header=T,sep=",",stringsAsFactors=F,fill=T)
  x.lumi <- lumiR(input.expressionfile,sep="\t")
  message ("Read x.lumi")
  lumi.N.Q <- lumiExpresso(x.lumi, normalize.param=list(method=norm))
  message ("Normalized Data")
  pdf("box-nonorm.pdf")
  boxplot(x.lumi)
  dev.off()
  pdf("hist-nonorm.pdf")
  hist(x.lumi)
  dev.off()
  pdf("box-norm.pdf")
  boxplot(lumi.N.Q)
  dev.off()
  pdf("hist-norm.pdf")
  hist(lumi.N.Q)
  dev.off()
  message ("plotted")
  pca<-prcomp(t(as.matrix(x.lumi)),scale=T)
  plot3d(pca$x[, 1:3], col=rainbow(nrow(phenod)), xlab="PC1", ylab = "PC2", zlab = "PC3", type = "s")
  rgl.postscript("PCA.pdf", fmt="pdf", drawText=TRUE)
  dataMatrix <- exprs(lumi.N.Q)
  presentCount <- detectionCall(x.lumi)
  if(filter.bool == TRUE) dataMatrix <- dataMatrix[presentCount > filter.int,]
  mData <- dataMatrix[, phenoline.start:phenoline.finish]
  mDataPheno<-phenod[phenoline.start:phenoline.finish, 1]
  design <- model.matrix(~ factor(mDataPheno, levels=c(phenoname.ctrl, phenoname.var)))
  fit <- lmFit(mData, design)
  fit <- eBayes(fit)
  tab <- topTable(fit, coef = 2, adjust = "fdr", n = dim(mData)[1])
  tab.sig <- tab[tab$adj.P.Val < pval & abs(tab$logFC) > abs(log2(fcval)),]
  tab.sig.in=tab.sig[rownames(tab.sig) %in% annotation$PROBE_ID,]
  mDataSig.in<-mData[rownames(mData) %in% rownames(tab.sig.in),]
  symbols=""
  for (i in 0:(nrow(mDataSig.in)) ) 
    symbols[i]=annotation$TargetID[annotation$PROBE_ID %in% rownames(mDataSig.in)[i]]
  pdf(heatmap.out)
  heatmap(mDataSig.in, labRow=symbols, col=topo.colors(16), cexCol=0.6)
  dev.off()
  symbols.table=""
  for (i in 0:(nrow(tab.sig.in)) ) 
    symbols.table[i]=annotation$TargetID[annotation$PROBE_ID %in% rownames(mDataSig.in)[i]]
  id.table=""
  for (i in 0:(nrow(tab.sig.in)) ) 
    id.table[i]=annotation$PROBE_ID[annotation$PROBE_ID %in% rownames(mDataSig.in)[i]]
  sequence.table=""
  for (i in 0:(nrow(tab.sig.in)) ) 
    sequence.table[i]=annotation$PROBE_SEQUENCE[annotation$PROBE_ID %in% rownames(mDataSig.in)[i]]
  definition.table=""
  for (i in 0:(nrow(tab.sig.in)) ) 
    definition.table[i]=annotation$DEFINITION[annotation$PROBE_ID %in% rownames(mDataSig.in)[i]]
  process.table=""
  for (i in 0:(nrow(tab.sig.in)) ) 
    process.table[i]=annotation$ONTOLOGY_PROCESS[annotation$PROBE_ID %in% rownames(mDataSig.in)[i]]
  function.table=""
  for (i in 0:(nrow(tab.sig.in)) ) 
    function.table[i]=annotation$ONTOLOGY_FUNCTION[annotation$PROBE_ID %in% rownames(mDataSig.in)[i]]
  component.table=""
  for (i in 0:(nrow(tab.sig.in)) ) 
    component.table[i]=annotation$ONTOLOGY_COMPONENT[annotation$PROBE_ID %in% rownames(mDataSig.in)[i]]
  tmp=data.frame(ID=id.table,Symbol=symbols.table,logFC=tab.sig.in$logFC,pVal=tab.sig.in$P.Value,adjPVal=tab.sig.in$adj.P.Val,Sequence=sequence.table,Definition=definition.table,Process=process.table,Function=function.table,Component=component.table)
  write.table(tmp, file=genes.out, quote=FALSE, row.names=FALSE, sep="\t")
  message ("Done")
}
