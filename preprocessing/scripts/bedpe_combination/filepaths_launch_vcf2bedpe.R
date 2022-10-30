library(data.table)

extract = function (folderpath){
  sink("Z:/wolu/testing_svaba/scripts/vcf2bedpe/samples.txt",append=TRUE)
  lof <- list.files(folderpath)
  log <- list.files("Z:/Shu/tcga_svaba_germline_marcin/")
  for (i in (1:length(lof))){
    split0 <- unlist(strsplit(lof[i],"[.]"))
    if (split0[5] == 'sv'){
      germline <- capture.output(cat(split0[1:2],"germline.sv.vcf",sep="."))
      if (germline %in% log){
        cat(paste0("/xchip/beroukhimlab/Shu/tcga_svaba_germline_marcin/",germline,"%/xchip/beroukhimlab/pcawgCalls/svaba/SN0246229/",lof[i]), sep = '\n', append = T)            
      }
      else{
        next
      }
    }
    else{
      next
    }    
  }
}

folder <- "Z:/siyun/data/insertions/pcawg/"
extract(folder)