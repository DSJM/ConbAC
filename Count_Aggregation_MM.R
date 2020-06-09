args <- commandArgs(TRUE)

source("SparseMatrix.R")

Count_Aggregation <- function(TrainDataDir, TestDataDir, OutputDir){
  "
  Count_Aggregation
  Function combines Train and Test Data to one csv, 
  also filters out features that are not present in both datasets.
  
  Parameters
  ----------
  Train/TestData : Train population count file path (.csv).
  OutputDir : Output directory defining the path of the exported files.
  "

  
  
  TrainData <- Read10x_New(TrainDataDir)
  
  TestData <- Read10x_New(TestDataDir)
  
  TrainData <- TrainData[intersect(rownames(TrainData), rownames(TestData)),]
  TestData <- TestData[intersect(rownames(TrainData), rownames(TestData)),]
  
  Data <- t(cbind(TrainData, TestData))
  
  Store10x_New(Data, paste0(OutputDir,'/data_combined'))
  
}


Count_Aggregation(args[1], args[2], args[3])
