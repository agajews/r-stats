rnaRandomTest=function(nclusters,ntests){
  predictionSum=0
  predictionCount=0
  for(i in 1:nclusters){
    prediction=rnaSeqSignature("crc_rnaseq.csv","phenoData.csv",",","upperquartile",0.4,longCluster,T,7,ntests,T,F)
    predictionSum=predictionSum+prediction
    predictionCount=predictionCount+1
  }
  predictionAve=predictionSum/predictionCount
  print(predictionAve)
  return(predictionAve)
}