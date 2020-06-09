
library(Matrix)

Read10x_New <- function(matrix_dir){
  barcode.path <- paste0(matrix_dir, "/barcodes.tsv")
  features.path <- paste0(matrix_dir, "/genes.tsv")
  matrix.path <- paste0(matrix_dir, "/matrix.mtx")
  mat <- readMM(file = matrix.path)
  feature.names = read.delim(features.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path, 
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = barcode.names$V1
  rownames(mat) = feature.names$V1
  
  return(mat)
}


Store10x_New <- function(Matrix,OutputDir){
  dir.create(OutputDir)
  writeMM(obj = Matrix, file=paste0(OutputDir,"/matrix.mtx"))
  write(x = rownames(Matrix),file=paste0(OutputDir,"/genes.tsv"))
  write(x = colnames(Matrix), file=paste0(OutputDir,"/barcodes.tsv"))
}
