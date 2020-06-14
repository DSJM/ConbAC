args <- commandArgs(TRUE)

Label_Aggregation <- function(TrainLabelP, TestLabelP, OutputDir){
  "
  CITE_Validation
  Function combines Train and Test set and returns train and test indices, 
  also filter out cell populations with less than 10 cells.
  It return a 'CV_folds.RData' file which then used as input to classifiers wrappers.
  
  Parameters
  ----------
  Train/TestLabels : Train population annotations file path (.csv).
  col_Index : column index (integer) defining which level of annotation to use,
  in case of multiple cell type annotations (default is 1)
  OutputDir : Output directory defining the path of the exported files.
  "

  
  TrainLabels <- as.matrix(read.csv(TrainLabelP, header = TRUE))
  TestLabels <- as.matrix(read.csv(paste0(TestLabelP, '/barcodes.tsv'), header = FALSE))
  TestLabels <- matrix("XXX", dim(TestLabels)[1], dim(TrainLabels)[2])
  TrainLabels <- cbind(TrainLabels, rep("Train", length(TrainLabels)))
  


  TestLabels <- cbind(TestLabels, rep("Test", length(TestLabels)))
  
  Labels <- rbind(TrainLabels, TestLabels)
  
  write.csv(x = Labels, row.names = FALSE,file = paste0(OutputDir, '/labfile.csv'))
  
}

Label_Aggregation(args[1], args[2], args[3])
