Seurat_to_MM <- function(SeuratObject, Label, Output_Label, OutputDir){
  
  DefaultAssay(object = SeuratObject) <- "RNA"
  
  write.csv(x = as.matrix(Label), file = Output_Label, append = FALSE, col.names = c("x"), row.names=FALSE)
  
  dir.create(OutputDir)
  obj = SeuratObject@assays$RNA@counts
  
  writeMM(obj = obj, file=paste0(OutputDir,"/matrix.mtx"))
  
  write(x = rownames(obj),file=paste0(OutputDir,"/genes.tsv"))
  write(x = colnames(obj), file=paste0(OutputDir,"/barcodes.tsv"))
}