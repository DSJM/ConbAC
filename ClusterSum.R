args <- commandArgs(TRUE)

ConfPath <- args[2]
OutputDir <- args[1]

clusterSum <- function(ConfPath){
  "
  Script to summarize all confusion matrices in one big for manual evaluation.
  "
  library(dplyr)
  library(abind)

  conf_files <- list.files(path=ConfPath, pattern="*.csv", full.names=TRUE, recursive=FALSE)
  
  conf_list <- sapply(conf_files, read.csv,header=TRUE, sep = ",", row.names = 1)
  names(conf_list) <- gsub(".*/","",names(conf_list))
  names(conf_list) <- gsub(".csv","",names(conf_list))
  
  cols <- unique(unlist(lapply(conf_list, names)))
  rows <- unique(unlist(lapply(conf_list, rownames)))
  
  conf_3D <- 
    matrix(NA, length(rows), length(cols)) %>% 
    `colnames<-`(cols) %>% 
    `rownames<-`(rows)
  
  conf_3D <- lapply(conf_list, function(df){
    conf_3D[rownames(df), names(df)] <- as.matrix(df)
    conf_3D
    }) %>% 
    list(along = 3) %>% 
    do.call(what = abind)
  
  
  return(conf_3D)
}


ConfPath <- gsub("/[^/]*$","",ConfPath)
results <- clusterSum(ConfPath)
saveRDS(results, file = file.path(OutputDir,"ClusterSummary.rds" ))
