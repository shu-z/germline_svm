library(data.table)
library(GenomicRanges)
library(gUtils)
library(parallel)
hg19_genes <- readRDS('/xchip/beroukhimlab/Alex/InProgress/inst/extdata/hg19_geneRanges.rds')
hg19_exons <- readRDS('/xchip/beroukhimlab/Alex/InProgress/inst/extdata/hg19_exonRanges.rds')

check_reformat <- function(i, df){
 
  row <- df[i,]
  if (row$chrom1 > row$chrom2) {
    return(cbind(chrom1=row$chrom2, start1=row$start2, end1=row$start2, chrom2=row$chrom1, 
                 start2=row$start1, end2=row$start1, name=row$name, score=row$score, 
                 strand1=row$strand2, strand2=row$strand1, row[,11:ncol(df)]))
  } else if (row$chrom1==row$chrom2 & row$start1>row$start2){
    return(cbind(chrom1=row$chrom2, start1=row$start2, end1=row$start2, chrom2=row$chrom1, 
                 start2=row$start1, end2=row$start1, name=row$name, score=row$score, 
                 strand1=row$strand2, strand2=row$strand1, row[, 11:ncol(df)]))
  } else {
    return(row)
  }
}

annot <- function(itter, bp) {
  sub <- bp[itter,]
  
  if(sub$chrom1 == sub$chrom2) {
    sub_gr <- GRanges(sub$chrom1, IRanges(as.numeric(as.character(sub$start1)), as.numeric(as.character(sub$start2))))
    
    genes_overlapped <- hg19_genes %&% sub_gr
    exons_overlapped <- hg19_exons %&% sub_gr
    
    genes_fract <- genes_overlapped %O% sub_gr
    
    if(any(genes_fract == 1)) {
      sub[, CN_annot := 1]
    } else {
      sub[, CN_annot := 0]
    }
    
    if(length(exons_overlapped) > 0) {
      sub[,exon_annot := 1]
    } else {
      sub[,exon_annot := 0]
    }
    
  } else {
    bp1_gr <- GRanges(sub$chrom1, IRanges(as.numeric(as.character(sub$start1)), width = 1))
    bp2_gr <- GRanges(sub$chrom2, IRanges(as.numeric(as.character(sub$start2)), width = 1))
    
    exons_bp1 <- hg19_exons %&% bp1_gr
    exons_bp2 <- hg19_exons %&% bp2_gr
    
    exons_impacted <- length(exons_bp1) + length(exons_bp2)
    
    if(exons_impacted > 0) {
      sub[, CN_annot := 0]
      sub[,exon_annot := 1]
    } else {
      sub[, CN_annot := 0]
      sub[,exon_annot := 0]
    }
  }
  sub[,svtype := ifelse(chrom1 == chrom2, ifelse(strand1 == strand2, ifelse(strand1 == '+', "h2hINV", "t2tINV"), ifelse(strand2 == "+", "DUP", "DEL")),"INTER")]
  return(sub)
}



filepath <- commandArgs(T)[1]
bedpe <- fread(paste0(filepath))

bedpe[chrom1 == 'X', chrom1 := '23']
bedpe[chrom2 == 'X', chrom2 := '23']
bedpe[chrom1 == "Y", chrom1 := '24']
bedpe[chrom2 == "Y", chrom2 := '24']
bedpe$chrom1 <- as.numeric(bedpe$chrom1)
bedpe$chrom2 <- as.numeric(bedpe$chrom2)
bedpe$start1 <- as.numeric(bedpe$start1)
bedpe$start2 <- as.numeric(bedpe$start2)


bedpe_formated <- rbindlist(lapply(1:nrow(bedpe), check_reformat, bedpe))
bedpe_formated$chrom1 <- as.character(bedpe_formated$chrom1)
bedpe_formated$chrom2 <- as.character(bedpe_formated$chrom2)

bedpe_formated[chrom1 == '23', chrom1 := 'X']
bedpe_formated[chrom2 == '23', chrom2 := 'X']
bedpe_formated[chrom1 == '24', chrom1 := 'Y']
bedpe_formated[chrom2 == '24', chrom2 := 'Y']

bedpe_annot <- rbindlist(mclapply(1:nrow(bedpe_formated), annot, bedpe_formated, mc.cores = 5))

sample_output <- paste0('/xchip/beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/data/4_annot/', unlist(strsplit(filepath,'3_fuzzy_match/'))[2], '.annot.bedpe')

write.table(bedpe_annot, sample_output, sep = '\t', row.names = F, col.names = T, quote = F)
