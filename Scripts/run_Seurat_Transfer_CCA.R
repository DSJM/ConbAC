args <- commandArgs(TRUE)

source("SparseMatrix.R")

run_Seurat_Transfer_CCA<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,GeneOrderPath = NULL,NumGenes = NULL){
  "
  run Seurat
  Wrapper script to run Seurat on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  Data <- Read10x_New(DataPath)
  Labels <- as.matrix(read.csv(LabelsPath))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                               Seurat                                      #
  #############################################################################
  library(Seurat)
  True_Labels_Seurat <- list()
  Pred_Labels_Seurat <- list()
  Total_Time_Seurat <- list()
  Data = t(as.matrix(Data))
  
  for (i in c(1:n_folds)){
    if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
      
      seur_Train <- CreateSeuratObject(counts = Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]])
      seur_Train@meta.data$Celltype = Labels[Train_Idx[[i]]]
      seur_Train <- NormalizeData(seur_Train)
      ## Limited to 5000 variable features to reduce calculation time, might be increased or removed
      seur_Train <- FindVariableFeatures(seur_Train, selection.method = "vst",nfeatures = NumGenes, verbose = FALSE)
      selected.genes <- rownames(seur_Train) #Take all Genes
      #selected.genes <- VariableFeatures(object = seur_Train) # Limit to variable genes
      seur_Train <- ScaleData(seur_Train, features = selected.genes)
      
      
      seur_Test <- CreateSeuratObject(counts = Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]])
      seur_Test@meta.data$Celltype = Labels[Test_Idx[[i]]]
      seur_Test <- NormalizeData(seur_Test)
      seur_Test <- ScaleData(seur_Test, features = selected.genes)
      
    }
    
    else{
      seur_Train <- CreateSeuratObject(counts = Data[,Train_Idx[[i]]])
      seur_Train@meta.data$Celltype = Labels[Train_Idx[[i]]]
      seur_Train <- NormalizeData(seur_Train)
      ## Limited to 5000 variable features to reduce calculation time, might be increased or removed
      seur_Train <- FindVariableFeatures(seur_Train, selection.method = "vst",nfeatures = 5000, verbose = FALSE)
      selected.genes <- rownames(seur_Train) #Take all Genes
      #selected.genes <- VariableFeatures(object = seur_Train) # Limit to variable genes
      seur_Train <- ScaleData(seur_Train, features = selected.genes)
      
      
      seur_Test <- CreateSeuratObject(counts = Data[,Test_Idx[[i]]])
      seur_Test@meta.data$Celltype = Labels[Test_Idx[[i]]]
      seur_Test <- NormalizeData(seur_Test)
      seur_Test <- ScaleData(seur_Test, features = selected.genes)
      
    }
    
    
    # Seurat Integration    
    start_time <- Sys.time()
    seur.anchors <- FindTransferAnchors(reference = seur_Train, query = seur_Test,reduction = "cca",features = selected.genes, dims = 1:30)
    predictions <- TransferData(anchorset = seur.anchors, refdata = seur_Train@meta.data$Celltype, weight.reduction = "cca",                                dims = 1:30)
    seur_Test <- AddMetaData(seur_Test, metadata = predictions)
    
    end_time <- Sys.time()
    
    Total_Time_Seurat[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    True_Labels_Seurat[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_Seurat[i] <- list(seur_Test@meta.data$predicted.id)
    
  }
  
  
  True_Labels_Seurat <- as.vector(unlist(True_Labels_Seurat))
  Pred_Labels_Seurat <- as.vector(unlist(Pred_Labels_Seurat))
  Total_Time_Seurat <- as.vector(unlist(Total_Time_Seurat))
  
  write.csv(True_Labels_Seurat,paste0(OutputDir,'/Seurat_CCA_true.csv'),row.names = FALSE)
  write.csv(Pred_Labels_Seurat,paste0(OutputDir,'/Seurat_CCA_pred.csv'),row.names = FALSE)
  write.csv(Total_Time_Seurat,paste0(OutputDir,'/Seurat_CCA_total_time.csv'),row.names = FALSE)
  }

if (args[6] == "0") {
  run_Seurat_Transfer_CCA(args[1], args[2], args[3], args[4])
} else {
  run_Seurat_Transfer_CCA(args[1], args[2], args[3], args[4], args[5], as.numeric(args[6]))
}