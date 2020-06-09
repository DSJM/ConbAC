withCallingHandlers({
  install.packages("devtools", repos="https://cloud.r-project.org/")
  install.packages("Seurat", repos="https://cloud.r-project.org/")
},
warning = function(w) stop(w))
