library(data.table)

extract = function (folderpath){
  sink("Z:/wolu/testing_svaba/scripts/fuzzymatch/filtered.txt",append=TRUE)
  lof <- list.files(folderpath)
  for (i in (1:length(lof))){
    cat(paste0("/xchip/beroukhimlab/wolu/testing_svaba/outputs/final_filter/filtered_only/",lof[i]),sep = "\n",append=T)
  }
}

folder <- "Z:/wolu/testing_svaba/outputs/final_filter/filtered_only/"
extract(folder)