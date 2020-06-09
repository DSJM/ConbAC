args <- commandArgs(TRUE)

ConfPath <- args[2]
OutputDir <- args[1]

clusterSum <- function(x){
  "
  Script to summarize all confusion matrices in one big for manual evaluation.
  "
  print(x)

}

clusterSum(ConfPath)

results <- "TEst"
saveRDS(results, file = file.path(OutputDir,"ClusterSummary.rds" ))
