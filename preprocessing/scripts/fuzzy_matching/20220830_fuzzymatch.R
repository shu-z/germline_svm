install.packages("/xchip/beroukhimlab/Ipsa/manta/scripts/InProgress.tar.gz", repos = NULL, type="source")

library(data.table)
library(InProgress)

run_it <- function(lof_pth) {
  sample.name <- unlist(strsplit(lof_pth, "2_filtered_data/"))[2]
  
  # output_path <- paste0('/xchip/beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/data/3_fuzzy_match/',sample.name,'_fuzzy.bedpe')
  output_path <- paste0('/xchip/beroukhimlab/siyun/germline_classifier/outputs/fuzzy_matching/',sample.name,'_fuzzy.bedpe')
  
  bedpe <- fread(lof_pth)
  bedpe[, chrom1 := gsub('chr','',chrom1)]
  bedpe[, chrom2 := gsub('chr','',chrom2)]
  
  fuzzy <- closest_germline(bp = bedpe, cores = 5, genome = 'hg19')
  line_sine <- closest_line_sine(bp = fuzzy, genome = 'hg19', cores = 5)
  
  write.table(line_sine, output_path, sep = '\t', row.names = F, col.names = T, quote = F)
}

filepath <- commandArgs(T)[1]
run_it(filepath)