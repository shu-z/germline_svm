library(data.table)

lof <- list.files('/Volumes/xchip_beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/data/1_bedpe_combined')
lof_pth <- paste0('/Volumes/xchip_beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/data/1_bedpe_combined/', lof) 


for(i in 1:length(lof_pth)) {
  cat(paste0(i,'\n'))
  tmp <- fread(paste0(lof_pth[i]))
  tmp[, EVDNC := gsub(".*?EVDNC=([A-Z]+).*", "\\1", INFO)]
  tmp[, NALT_SR := unlist(strsplit(TUMOR, ":"))[3], by = 'TUMOR']
  tmp[, NALT := unlist(strsplit(TUMOR, ":"))[1], by = 'TUMOR']
  
  tmp_mq <- tmp[MAPQ_1 == 60 & MAPQ_2 == 60]
  tmp_mq_evd <- tmp_mq[!(EVDNC == 'DSCRD')]
  tmp_fin <- tmp_mq_evd[NALT_SR > 1]
  tmp_span <- tmp_fin[SPAN > 49 | SPAN== -1]
  write.table(tmp_span, paste0('/Volumes/xchip_beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/data/2_filtered_data/', lof[i],'.filtered'),
              sep = '\t', row.names = F, col.names = T, quote = F)
}
