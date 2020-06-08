args <- commandArgs(TRUE)

Annot_Fold <- function(LabelPath, col_Index = 1,OutputDir, minRefCells = 10){
  "
  Annot_Fold
  Function takes combined label set and returns 
  1. train and test indices, 
  2. filter out ref cell populations with less than 10 cells.
  3. a 'CV_folds.RData' file which then used as input to classifiers wrappers.
  
  Parameters
  ----------
  LabelPath : Combined population annotations file path (.csv).
  col_Index : column index (integer) defining which level of annotation to use,
  in case of multiple cell type annotations (default is 1)
  OutputDir : Output directory defining the path of the exported file.
  "
  
  Labels <- read.csv(LabelPath)
  
  ## Removing cells that are not sufficient present in Reference (or Test cluster)
  Removed_classes <- !(table(Labels[Labels[, -1] == "Train",][,col_Index]) > minRefCells | table(Labels[Labels[, -1] == "Test",][,col_Index]) >= 1)
  Cells_to_Keep <- !(is.element(Labels[,col_Index],names(Removed_classes)[Removed_classes]))
  Labels <- Labels[Cells_to_Keep,]

  
  ## Retrieving CV_folds_variables: 
  n_folds = 1
  Train_Idx <- list(which(Labels[, -1] == "Train"))
  Test_Idx <- list(which(Labels[, -1] == "Test"))
  
   
  save(n_folds,Train_Idx,Test_Idx,col_Index,Cells_to_Keep,file = paste0(OutputDir, '/CV_folds.RData'))
}

Annot_Fold(args[1], as.numeric(args[2]), args[3])