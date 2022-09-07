lof <- list.files('/Volumes/xchip_beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/data/2_filtered_data')
lof_pth <- paste0('/xchip/beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/data/2_filtered_data/', lof)

sink('/Volumes/xchip_beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/scripts/fuzzy_matching/20220418_TApaths.sh')
for(i in 1:length(lof_pth)) {
  cat(paste0('Rscript ', 
             '/xchip/beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/scripts/fuzzy_matching/20220418_fuzzymatch_snowman.R ',
             lof_pth[i]), sep = '\n', append = T)
}
