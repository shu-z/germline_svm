lof <- list.files('/Volumes/xchip_beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/data/3_fuzzy_match')
lof_pth <- paste0('/xchip/beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/data/3_fuzzy_match/', lof)

sink('/Volumes/xchip_beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/scripts/annot/20220419_annot_TA.sh')
for(i in 1:length(lof_pth)) {
  cat(paste0("Rscript ",
             "/xchip/beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/scripts/annot/20220419_annot.R ",
             lof_pth[i]), append = T, sep = '\n')
}
closeAllConnections()
