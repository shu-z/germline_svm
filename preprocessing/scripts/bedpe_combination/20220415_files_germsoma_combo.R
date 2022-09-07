require(data.table)
lof_pths <- list.files('/Volumes/xchip_beroukhimlab/Alex/CCLE_PCAWG/Germline_Correctiion/Test_Set')

ids_pres <- sapply(lof_pths, function(x) {unlist(strsplit(x, "_"))[1]})

lof_germline <- list.files('/Volumes/xchip_beroukhimlab/Simona/SV_Homeology/Data/pdcData/actualData')
lof_germline <- lof_germline[(grepl("germline", lof_germline))]
lof_germline <- lof_germline[(grepl("broad-snowman", lof_germline))]

lof_somatic <- list.files('/Volumes/xchip_beroukhimlab/Simona/SV_Homeology/Data/pdcData/actualData')
lof_somatic <- lof_somatic[(grepl("somatic", lof_somatic))]
lof_somatic <- lof_somatic[(grepl("broad-snowman", lof_somatic))]
lof_somatic <- lof_somatic[!(grepl("vcf.gz", lof_somatic))]

lof_germ_need_pth <- paste0("/xchip/beroukhimlab/Simona/SV_Homeology/Data/pdcData/actualData/", lof_germline)
lof_soma_need_pth <- paste0("/xchip/beroukhimlab/Simona/SV_Homeology/Data/pdcData/actualData/", lof_somatic)
c <- as.data.table(cbind(lof_germ_need_pth, lof_soma_need_pth))
cc <- paste0(c$lof_germ_need_pth, "%", c$lof_soma_need_pth)
sink('/Volumes/xchip_beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/scripts/bedpe_combination/20220415_germsoma_sublist.sh')

for(i in 1:length(cc)) {

  cat(paste0("Rscript", " ", 
             "/xchip/beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/scripts/bedpe_combination/20220415_vcf2bedpe_andcombine.R", " ",
             cc[i]), 
      file = '/Volumes/xchip_beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/scripts/bedpe_combination/20220415_germsoma_sublist.sh',
      sep = '\n', append = TRUE)
  
}
closeAllConnections()
