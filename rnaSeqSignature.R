rnaSeqSignature = function(input.table,pheno.table,sep,norm,Delta,v.cluster,balance,test.set,tests,randomTest,t.stats){
  #source("http://bioconductor.org/biocLite.R")
  #biocLite("cqn")
  #biocLite("edgeR")
  require(edgeR)
  require(baySeq)
  require(pheatmap)
  require(DESeq2)
  require(EnsDb.Mmusculus.v79)
  require(annotate)
  require(tools)
  require(pamr)
  rawCountTable=read.table(input.table,header=T,stringsAsFactors=F,row.names=1,sep=sep)
  rawCountTable=rawCountTable[complete.cases(rawCountTable),]
  rawCountTable=as.matrix(rawCountTable)
  if(randomTest==T){
    test.set.three=test.set+3
    randClusterSize=sample(test.set.three:(ncol(rawCountTable)-test.set.three), 1)
    v.cluster=colnames(rawCountTable)[sample(1:ncol(rawCountTable),randClusterSize,replace=FALSE)]
    print(v.cluster)
    #return(v.cluster)
  }
  longCount=rawCountTable[,colnames(rawCountTable) %in% v.cluster]
  otherCount=rawCountTable[,!(colnames(rawCountTable) %in% v.cluster)]
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
  phenoTable=read.table(pheno.table,header=T,stringsAsFactors=F,sep="\t")
  pseudoCount = log2(rawCountTable + 1)
  isCount = rowSums(pseudoCount != 0)
  filteredCount = pseudoCount[isCount > 0.8*ncol(rawCountTable),]
  filteredDGE=DGEList(counts=filteredCount,group=factor(names(filteredCount)))
  print(paste(nrow(filteredDGE),"genes after filtering"))
  normCount=calcNormFactors(filteredDGE,method=norm,p=0.75)
  normClass = vector()
  for(i in 1:ncol(normCount)){
    if(colnames(normCount)[i] %in% colnames(longCount)){
      normClass[i]="long"
    }
    else if(colnames(normCount)[i] %in% colnames(otherCount)){
      normClass[i]="other"
    }
  }
  #return(normCount)
  #classCount=data.frame(classCount)
  #classCount=classCount[,order(classCount$Class)]
  normCount=as.matrix(normCount)
  #print(paste("Norm Count Dim:",dim(normCount)))
  Class=normClass
  classCount=rbind(colnames(normCount),Class)
  #return(classCount)
  print("Merged")
  #rankedList=getRLs(as.matrix(otherCount),as.matrix(longCount))
  #return(rankedList)
  #return(classCount
  #classCount=t(classCount)
  #normCount=t(normCount)
  sdevs <- apply(normCount, 1, sd)
  eset <- normCount[sdevs > median(sdevs), ]
  #t.stats <- apply(classCount, 1, function(e) {
    #t.test(e ~ classCount[,nrow(classCount)])$statistic
  #})
  #ranks <- order(abs(t.stats), decreasing = TRUE)
  #selected.exprs <- classCount[ranks[1:50], ]
  #return(selected.exprs)
  #dds=DESeq(normCount)
  #return(eset)
  #return(known.labels)
  if(t.stats==T){
    print("Running T.Stats")
    t.stats <- apply(eset, 1, function(e){t.test(e ~ classCount[2,])$statistic})
    ranks <- order(abs(t.stats), decreasing = TRUE)
    selected.exprs <- eset[ranks[1:50], ]
  }
  print("Getting Symbols")
  
  
  #symbols.norm=c()
  #rownames=rownames(normCount)
  ##return(known.patients)
  #for(i in 1:length(rownames)){
  #  if(rownames[i]!="ENSG00000280804" && rownames[i]!="ENSG00000280804"){
  #    symbol=gsub("\\..*","", rownames[i])
  #    symbolmatrix=select(org.Mm.eg.db,keys=symbol,columns=c("ENSEMBL","SYMBOL","GENENAME"),keytype="ENSEMBL")
  #    symbols.norm[i] <- symbolmatrix$ENSEMBL
  #    if(symbol=='' || is.null(symbol)){
  #      symbol='NA'
  #      symbols.norm[i] <- symbol[1,1]
  #    }
  #  }
  #  else{
  #    symbol='NA'
  #    symbols.norm[i] <- symbol$ENSEMBL
  #  }
  #}
  #normSymbols=cbind(normCount,symbols.norm)
  #return(normSymbols)
  
  #symbols=c()
  #rownames=rownames(selected.exprs)
  #features.selected=rownames(selected.exprs)
  #for(i in 1:length(rownames)){
  #  symbols[i] <- getGene(gsub("\\..*","", rownames[i]))$symbol
  #}
  #rownames(selected.exprs)=symbols
  
  #symbols.known=c()
  #features.known=rownames(known.patients)
  #rownames=rownames(known.patients)
  #return(known.patients)
  #for(i in 1:length(rownames)){
  #  if(rownames[i]!="ENSG00000280804"){
  #    symbol=getGene(gsub("\\..*","", rownames[i]))$symbol
  #    if(symbol=='' || is.null(symbol)){
  #      symbol='NA'
  #    }
  #  }
  #  else{
  #    symbol='NA'
  #  }
  #  symbols.known[i] <- symbol
  #}
  #rownames(known.patients)=symbols.known
  
  #symbols.new=c()
  #features.new=rownames(new.patients)
  #rownames=rownames(new.patients)
  #for(i in 1:length(rownames)){
  #  if(rownames[i]!="ENSG00000280804"){
  #    symbol=getGene(gsub("\\..*","", rownames[i]))$symbol
  #    if(symbol=='' || is.null(symbol)){
  #      symbol='NA'
  #    }
  #  }
  #  else{
  #    symbol='NA'
  #  }
  #  symbols.new[i] <- symbol
  #}
  #rownames(new.patients)=symbols.new
  testSignature = function(){
    hide=eset[,sample(1:ncol(eset),test.set,replace=FALSE)]
    known.patients <- eset[, !(colnames(eset) %in% colnames(hide))]
    new.patients <- eset[, colnames(eset) %in% colnames(hide)]
    print("Selected Patients")
    known.labels <- classCount[2,classCount[1,] %in% colnames(known.patients)]
    new.labels <- classCount[2,classCount[1,] %in% colnames(new.patients)]
    features.known=rownames(known.patients)
    #rownames(known.patients)=normSymbols[rownames(normSymbols) %in% rownames(known.patients),ncol(normSymbols)]
    dat=known.patients
    gN=rownames(known.patients)
    gI=features.known
    sI=colnames(known.patients)
    print("Training Model...")
    train.dat <- list(x = dat, y = known.labels, geneid = gI, sampleid = sI)
    model <- pamr.train(train.dat, n.threshold = 100)
    if(balance==F){
      model.cv <- pamr.cv(model, train.dat, nfold = 10)
      pamr.plotcv(model.cv)
    }
    mean.class.error <- function(model, threshold, good.error.weight = 0.5) {
      conf.matrix <- pamr.confusion(model, threshold, FALSE)
      good.error <- conf.matrix["good", "poor"]/sum(conf.matrix["good",])
      poor.error <- conf.matrix["poor", "good"]/sum(conf.matrix["poor", ]) 
      good.error.weight * good.error + (1 - good.error.weight) * poor.error}
    plot.mean.class.error <- function(model, good.error.weight = 0.5) {
      weighted.errors <- sapply(model$threshold, function(t) {
        mean.class.error(model, t, good.error.weight)
      })
      plot(model$threshold, weighted.errors, type = "l", col = 4, ylim = c(-0.1, 1.1), ylab = "Mean class error", xlab = "Threshold") 
      axis(3, at = model$threshold, lab = paste(model$size), srt = 90, adj = 0)
      mtext("Number of genes", 3, 3)
      grid()
    }
    #return(model.cv)
    pamr.predict = pamr.predict(model, new.patients, Delta)
    pamr.predict.names = rownames(pamr.predict(model, new.patients, Delta,type="posterior"))
    prediction=cbind(pamr.predict.names,as.vector(pamr.predict),new.labels)
    colnames(prediction)=c("Sample","Prediction","Class")
    print("\n")
    print(prediction)
    return(prediction)
  }
  testsCorrect=0
  testsWrong=0
  totalTests=0
  for(i in 1:tests){
    prediction=testSignature()
    for(i in 1:test.set){
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
  #pamr.listgenes(model,train.dat,Delta)
  print(predictability)
  return(predictability)
}