library(data.table)

lof <- list.files('/xchip/beroukhimlab/wolu/testing_svaba/outputs/final_bedpe/final_bedpe_only') #folder containing the bedpe files
lof_pth <- paste0('/xchip/beroukhimlab/wolu/testing_svaba/outputs/final_bedpe/final_bedpe_only/', lof) 


for(i in 1:length(lof_pth)) {
  cat(paste0(i,'\n'))
  tmp <- fread(paste0(lof_pth[i]))
  tmp[, EVDNC := gsub(".*?EVDNC=([A-Z]+).*", "\\1", INFO)]
  tmp[, NALT_SR := unlist(strsplit(TUMOR, ":"))[3], by = 'TUMOR']
  tmp[, NALT := unlist(strsplit(TUMOR, ":"))[1], by = 'TUMOR']
  
  #tmp_mq <- tmp[MAPQ_1 == 60 & MAPQ_2 == 60] #takes only reads that have a mapq score that are 60 or higher since MAPQ scores are capped at 60
  tmp_mq <- tmp[MAPQ_1 == 60 | MAPQ_2 == 60] #takes only reads that have a mapq score that are 60 or higher since MAPQ scores are capped at 60
  tmp_mq_evd <- tmp_mq[!(EVDNC == 'DSCRD')] #removes reads that are discordant
  tmp_fin <- tmp_mq_evd[NALT_SR > 1] #the number of reads covering the site/depth of coverage must be greater than 1
  tmp_span <- tmp_fin[SPAN > 49 | SPAN== -1] #SPAN = -1 refers to a translocation. Span shorter than 50bp are considered simple indels
  write.table(tmp_span, paste0('/xchip/beroukhimlab/wolu/testing_svaba/outputs/final_filter/filtered_only/', lof[i],'.filtered'),
              sep = '\t', row.names = F, col.names = T, quote = F)
}
