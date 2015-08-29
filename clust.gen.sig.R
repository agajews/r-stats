clust.gene.sig = function(input.table,input.sep="\t",out.sep="\t",v.clust,cluster.label="Selected",other.label="Other",norm="upperquartile",threshold=0.9,balance=F,test.size=7,ntests=50,random=F,ngenes=50,return=F,preprocess=T,fold.change=2,max.adj.p=0.05,min.expr=0.8,symbols=T){
  require(edgeR)
  require(tools)
  require(pamr)
  require(mygene)
  print("Reading Inputs...")
  rawCountTable=read.table(input.table,header=T,stringsAsFactors=F,row.names=1,sep=input.sep,fill=T)
  rawCountTable=rawCountTable[complete.cases(rawCountTable),]
  rawCountTable=as.matrix(rawCountTable)
  if(ncol(rawCountTable)==0 || ncol(rawCountTable)==1){
    print("Incorrect input.sep value")
    return()
  }
  if(v.clust=='NA'){
    print("Error: Please define a cluster or enable random cluster selection")
    return()
  }
  if(random==T){
    test.size.three=test.size+3
    randClusterSize=sample(test.size.three:(ncol(rawCountTable)-test.size.three), 1)
    v.clust=colnames(rawCountTable)[sample(1:ncol(rawCountTable),randClusterSize,replace=FALSE)]
  }
  longCount=rawCountTable[,colnames(rawCountTable) %in% v.clust]
  otherCount=rawCountTable[,!(colnames(rawCountTable) %in% v.clust)]
  if(balance==TRUE && ncol(longCount)<ncol(otherCount)){
    otherSel=otherCount[,sample(1:ncol(otherCount),ncol(longCount),replace=FALSE)]
  }
  else{
    otherSel=otherCount[,sample(1:ncol(otherCount),ncol(otherCount),replace=FALSE)]
  }
  otherSel=as.matrix(otherSel)
  longCount=as.matrix(longCount)
  print(colnames(longCount))
  rawCountTable=data.frame(longCount,otherSel)
  if(preprocess==T){
    pseudoCount = log2(rawCountTable + 1)
    isCount = rowSums(pseudoCount != 0)
    filteredCount = pseudoCount[isCount > min.expr*ncol(rawCountTable),]
    filteredDGE=DGEList(counts=filteredCount,group=factor(names(filteredCount)))
    print(paste(nrow(filteredDGE),"genes after filtering"))
    normCount=calcNormFactors(filteredDGE,method=norm,p=0.75)
  }
  else{
    normCount=rawCountTable
  }
  #return(normCount)
  normClass = vector()
  for(i in 1:ncol(normCount)){
    if(colnames(normCount)[i] %in% colnames(longCount)){
      normClass[i]=cluster.label
    }
    else if(colnames(normCount)[i] %in% colnames(otherCount)){
      normClass[i]=other.label
    }
  }
  normCount=as.matrix(normCount)
  Class=normClass
  classCount=rbind(colnames(normCount),Class)
  sdevs <- apply(normCount, 1, sd)
  eset <- normCount[sdevs > median(sdevs), ]
  t.stats=F
  if(t.stats==T){
    print("Running T.Stats...")
    t.stats <- apply(eset, 1, function(e){t.test(e ~ classCount[2,])$statistic})
    #return(t.stats)
    ranks <- order(abs(t.stats), decreasing = TRUE)
    eset <- eset[ranks[1:ngenes], ]
  }
  if(t.stats==F){
    print("Selecting Genes...")
    mDataPheno=vector()
    for(i in 1:ncol(normCount)){
      mDataPheno[i]=classCount[2,classCount[1,] %in% colnames(normCount)[i]]
    }
    mData=normCount
    design <- model.matrix(~ factor(mDataPheno, levels=c(other.label, cluster.label)))
    colnames(design) <- c(other.label, paste(other.label,".vs.",cluster.label,sep=""))
    fit=lmFit(mData,design)
    fit=eBayes(fit)
    tab=topTable(fit,coef=2,adjust="fdr",n=dim(mData)[1])
    tab.sig=tab[1,]
    rowCount=1
    for(i in 1:nrow(tab)){
      if(tab$adj.P.Val[i]>max.adj.p){
        if(2^abs(tab$logFC[i])>fold.change){
          tab.sig[rowCount,]=tab[i,]
          rowCount=rowCount+1
        }
      }
    }
    Fold.Change=vector()
    for(i in 1:nrow(tab.sig)){
      if(tab.sig$logFC[i]>0){
        Fold.Change[i]=2^(abs(tab.sig$logFC[i]))
      }
      else{
        Fold.Change[i]=(-1)*(2^(abs(tab.sig$logFC[i])))
      }
    }
    tab.sig=cbind(tab.sig,Fold.Change)
    if(nrow(tab.sig)>2){
      print(paste(nrow(tab.sig),"Differentially Expressed Genes"))
      eset=eset[rownames(eset) %in% as.vector(rownames(tab.sig)),]
    }
  }
  if(nrow(tab.sig)<3){
    print("Insufficient Significant Genes Detected With Given 'fold.change' & 'max.adj.p' Parameters")
    print(paste("Running T.Stats To Compute Top",ngenes,"Differentially Expressed Genes Instead..."))
    t.stats <- apply(eset, 1, function(e){t.test(e ~ classCount[2,])$statistic})
    ranks <- order(abs(t.stats), decreasing = TRUE)
    eset <- eset[ranks[1:ngenes], ]
  }
  if(symbols==T){
    print("Getting Symbols...")
    symbols=c()
    rownames=rownames(eset)
    features.selected=rownames(eset)
    for(i in 1:length(rownames)){
      symbol=getGene(gsub("\\..*","", rownames[i]))$symbol
      if(is.null(symbol)){symbol="NA"}
      symbols[i]=symbol
    }
    rownames(eset)=symbols
    if(nrow(tab.sig)>2){
      Gene_Symbols=symbols
      sig.genes=cbind(Gene_Symbols,tab.sig)
      write.table(sig.genes,file="clust.gene.sig.csv",row.names=F,sep=out.sep)
    }
  }
  testSignature = function(){
    print("Selecting Patients...")
    hide=eset[,sample(1:ncol(eset),test.size,replace=FALSE)]
    known.patients <- eset[, !(colnames(eset) %in% colnames(hide))]
    new.patients <- eset[, colnames(eset) %in% colnames(hide)]
    known.labels <- classCount[2,classCount[1,] %in% colnames(known.patients)]
    new.labels <- classCount[2,classCount[1,] %in% colnames(new.patients)]
    features.known=rownames(known.patients)
    dat=known.patients
    gN=rownames(known.patients)
    gI=features.known
    sI=colnames(known.patients)
    print("Training Model...")
    train.dat <- list(x = dat, y = known.labels, genenames = gN, geneid = gI, sampleid = sI)
    model <- pamr.train(train.dat, n.threshold = 100)
    if(balance==F){
      model.cv <- pamr.cv(model, train.dat, nfold = 10)
      pamr.plotcv(model.cv)
    }
    pamr.predict = pamr.predict(model, new.patients, threshold)
    pamr.predict.names = rownames(pamr.predict(model, new.patients, threshold,type="posterior"))
    prediction=cbind(pamr.predict.names,as.vector(pamr.predict),new.labels)
    colnames(prediction)=c("Sample","Prediction","Class")
    print("\n")
    #genes.list=pamr.listgenes(model,train.dat,threshold,genenames=T)
    print(prediction)
    return(prediction)
  }
  testsCorrect=0
  testsWrong=0
  totalTests=0
  for(i in 1:ntests){
    prediction=testSignature()
    for(i in 1:test.size){
      if(prediction[i,2]==prediction[i,3]){
        testsCorrect=testsCorrect+1
      }
      else{
        testsWrong=testsWrong+1
      }
      totalTests=totalTests+1
    }
  }
  predictability=testsCorrect/totalTests
  if(nrow(tab.sig)>2){
    numDiffGenes=nrow(eset)
    print(paste(numDiffGenes,"Differentially Expressed Genes"))
  }
  print(paste("Predictability of Cluster:",predictability))
  if(return==T){return(predictability)}
}

clust.gene.sig.rand=function(input.table,input.sep="\t",threshold=0.9,norm="upperquartile",test.size=7,nclusters=5,ntests=5,ngenes=50,return=F,fold.change=2,max.adj.p=0.05){
  predictionSum=0
  predictionCount=0
  for(i in 1:nclusters){
    prediction=clust.gene.sig(input.table=input.table,symbols=F,input.sep=input.sep,norm=norm,threshold=threshold,v.clust=c(""),balance=T,test.size=test.size,ntests=ntests,fold.change=fold.change,max.adj.p=max.adj.p,random=T,ngenes=ngenes,return=T)
    predictionSum=predictionSum+prediction
    predictionCount=predictionCount+1
  }
  predictionAve=predictionSum/predictionCount
  print(paste("Average Predictability of Random Clusters:",predictionAve))
  if(return==T){return(predictionAve)}
}

unsupervised.clust.gen = function(input.table,input.sep="\t",pheno.table,pheno.sep="\t",table.out="unsupervised.clust.gen.csv",norm="upperquartile",clusters=15,recursive=T,type="hclust",return=F){
  require(edgeR)
  require(baySeq)
  require(pheatmap)
  rawCountTable=read.table(input.table,header=T,stringsAsFactors=F,row.names=1,sep=input.sep)
  phenoTable=read.table(pheno.table,header=T,stringsAsFactors=F,sep=pheno.sep)
  phenoTable=phenoTable[sort(rownames(phenoTable)),]
  pseudoCount = log2(rawCountTable + 1)
  isCount = rowSums(pseudoCount != 0)
  filteredCount = pseudoCount[isCount > 0.8*ncol(rawCountTable),]
  filteredDGE=DGEList(counts=filteredCount,group=factor(names(rawCountTable)))
  print(paste(nrow(filteredDGE),"genes after filtering"))
  normCount=calcNormFactors(filteredDGE,method=norm,p=0.75)
  normMatrix = t(as.matrix(normCount))
  rnaClusters=rownames(normMatrix)
  if(recursive==T){
    for (i in 2:clusters){
      if(type=="kmeans"){
        name=paste("unsupervised.kmeans.set",i,".png",sep="")
        png(name)
        rnaHeatmap = pheatmap(normMatrix, kmeans_k = i, show_colnames=F, clustering_method = Ward.D2)
        dev.off()
        rnaClusters = cbind(rnaClusters, clusters = rnaHeatmap$kmeans$cluster)
      }
      else if(type=="hclust"){
        name=paste("unsupervised.hclust.tree",".png",sep="")
        d = dist(normMatrix,method="euclidean")
        hclust=hclust(d,"ward.D2")
        png(name)
        plot(hclust)
        dev.off()
        hcut=cutree(hclust,k=i)
        rnaClusters = cbind(rnaClusters, clusters = hcut)
      }
      else{
        print("Unknown Type")
      }
      print(paste("Finished Set",i))
    }
    colnames(rnaClusters)[2:ncol(rnaClusters)] = paste(2:clusters,"_Clusters",sep="")
    colnames(rnaClusters)[1]="Sample"
    rnaClusters=rnaClusters[sort(rownames(rnaClusters)),]
    writeCols=rnaClusters[,i:clusters]
  }
  else{
    if(type=="kmeans"){
      name=paste("unsupervised.kmeans.set",clusters,".png",sep="")
      png(name)
      rnaHeatmap = pheatmap(normMatrix, kmeans_k = clusters, show_colnames=F, clustering_method = Ward.D2)
      dev.off()
      rnaClusters = cbind(rnaClusters, clusters = rnaHeatmap$kmeans$cluster)
    }
    else if(type=="hclust"){
      name=paste("unsupervised.hclust.tree",".png",sep="")
      d = dist(normMatrix,method="euclidean")
      hclust=hclust(d,"ward.D2")
      png(name)
      plot(hclust)
      dev.off()
      hcut=cutree(hclust,k=clusters)
      rnaClusters = cbind(rnaClusters, clusters = hcut)
    }
    else{
      print("Unknown Type")
    }
    colnames(rnaClusters)[2] = paste(clusters,"_Clusters",sep="")
    colnames(rnaClusters)[1]="Sample"
    rnaClusters=rnaClusters[sort(rownames(rnaClusters)),]
    writeCols=rnaClusters[,2]
    print("Finished")
  }
  phenoTableOut=cbind(phenoTable,rnaClusters)
  write.table(phenoTableOut, file=table.out, quote=F, row.names=F, col.names=T, sep="\t")
  if(return==T){return(phenoTable)}
}