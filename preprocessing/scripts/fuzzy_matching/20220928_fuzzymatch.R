library(gUtils)
library(data.table)
library(GenomicRanges)
library(parallel)
#library(InProgress)   ##commenting out this line and all package-derived files; adding the individual functions and files here
#This script does the vcf to bedpe conversion, fuzzy filtering, and line/sine calculation.
##adding closest_germline function from Alex's fuzzy_filter.R

chrom1=chrom2=str_dist=end_dist=gnomad_dist=start=end=gnomad_germline=NULL
#' gNOMAD dataset to fuzzy filter
#'
#' gNOMAd v2.1 control sites lifted over to hg38
#' @name gnomad_germline_hg38all
#' @docType data
#' @keywords data internal
#' @format \code{data.table}
#gnomad_germline_hg38all = fread(system.file('extdata', 'gnomad_germline_hg38all.txt', package = 'InProgress'))
gnomad_germline_hg38all = fread("/xchip/beroukhimlab/wolu/testing_svaba/extdata/gnomad_germline_hg38all.txt") #changing these to specific file paths so that we get rid of the 'InProgress' reference

#' gNOMAD dataset to fuzzy filter
#'
#' gNOMAd v2.1 control sites in hg19
#' @name gnomad_germline_hg19all
#' @docType data
#' @keywords data internal
#' @format \code{data.table}
#gnomad_germline_hg19all = fread(system.file('extdata', 'gnomad_germline_hg19all.txt', package = 'InProgress'))
gnomad_germline_hg19all = fread("/xchip/beroukhimlab/wolu/testing_svaba/extdata/gnomad_germline_hg19all.txt")


#' @name fuzzy_filter_germline
#' @title Distance to closest germline annotator
#' @param itter passed from mclapply to iterate
#' @param bed Bedpe returned from annotate_sv function
#' @param g genome passed
#' @return SV data table with columns added indicating germline or somatic, germline is defined as <=1kbp away from agnostic perfect match in reference
#' @description 
#' 
#' Determines if each SV should be considered germline by hard filtering. Used with wrapper function.
#' 
#' @import data.table 
#' @keywords internal
fuzzy_filter_germline = function(itter = NULL, bed = NULL, g = NULL) {
  sub <- bed[itter,]
  ## reorder for filtering
  if(sub$chrom1 > sub$chrom2) {
    sub_ord <- cbind(chrom1=sub$chrom2, start1=sub$start2, end1=sub$end2, chrom2=sub$chrom1, start2=sub$start1, end2=sub$end1, sub[,7:ncol(bed)])
  } else if (sub$chrom1 == sub$chrom2 & sub$start1 > sub$start2) {
    sub_ord <- cbind(chrom1=sub$chrom2, start1=sub$start2, end1=sub$end2, chrom2=sub$chrom1, start2=sub$start1, end2=sub$end1, sub[,7:ncol(bed)])
  } else {
    sub_ord <- sub
  }
  ### change to integers to match reference germline
  sub_ord[chrom1 == "X", chrom1 := 23]
  sub_ord[chrom2 == "X", chrom2 := 23]
  sub_ord[chrom1 == "Y", chrom1 := 24]
  sub_ord[chrom2 == "Y", chrom2 := 24]
  ### subset reference to matching chromosome
  ref_sub <- g[chrom1 == sub_ord$chrom1 & chrom2 == sub_ord$chrom2]
  ### calculate distances
  ref_sub[,str_dist := abs(start - as.numeric(as.character(sub_ord$start1)))]
  ref_sub[,end_dist := abs(end - as.numeric(as.character(sub_ord$start2)))]
  ref_sub[,gnomad_dist := (str_dist + end_dist)]
  ### choose closest match
  ref_min <- ref_sub[which.min(ref_sub$gnomad_dist)]
  if(nrow(ref_min) == 0) {
    sub <- cbind(sub, gnomad_dist = "")
    return(sub)
  } else {
    sub <- cbind(sub, gnomad_dist = ref_min$gnomad_dist)
    return(sub)
  }
}

#' @name closest_germline
#' @title Determines distance to nearest germline event
#' @param bp \href{https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format}{Bedpe} from \link[InProgress]{svaba_vcf2bedpe} or \link[InProgress]{manta_vcf2bedpe}
#' @param cores Number of cores to run on, default is 1
#' @param genome run under hg19 or hg38
#' @return \href{https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format}{Bedpe} with a column added for distance to nearest germline event
#' @description 
#' 
#' Uses \href{https://gnomad.broadinstitute.org/downloads#v2-structural-variants}{gnomAD} to annotate the nearest germline event to each structural variant.
#' For more information read \href{https://www.nature.com/articles/s41586-020-2287-8}{gnomAD blog}. Reference is in hg38.
#' 
#' @import data.table
#' @importFrom parallel mclapply 
#' @export
closest_germline = function(bp = NULL, cores = 1, genome = NULL) {
  if(is.null(bp)) {
    stop('NULL input')
  }
  
  if(as.character(genome) == 'hg19') {
    gnomad_germline = gnomad_germline_hg19all
  } else if (as.character(genome) == 'hg38') {
    gnomad_germline = gnomad_germline_hg38all
  } else {
    stop('Please state hg19 or hg38 as genome')
  }
  
  cat("Comparing against known germline...")
  annotated_bedpe <- rbindlist(mclapply(1:nrow(bp), fuzzy_filter_germline, bp, g = gnomad_germline, mc.cores = cores))
  cat("done.\n")
  return(annotated_bedpe)
}

##adding closest_line_sine function from Alex's line_sine_annotation.R script
str1_dist_l=str1_dist_s=str2_dist_l=str2_dist_s=line_dist=sine_dist=LINE_dt=SINE_dt=NULL
#' Line repeatmaskers
#'
#' repeat masker line elements in hg38
#' @name LINE_dt_hg38
#' @docType data
#' @keywords data internal
#' @format \code{data.table}
#LINE_dt_hg38 = readRDS(system.file('extdata', 'repeatmasker_hg38_LINE.bed', package = 'InProgress'))
LINE_dt_hg38 = readRDS("/xchip/beroukhimlab/wolu/testing_svaba/extdata/repeatmasker_hg38_LINE.bed") 


#' Sine repeatmaskers
#'
#' repeat masker sine elements in hg38
#' @name SINE_dt_hg38
#' @docType data
#' @keywords data internal
#' @format \code{data.table}
#SINE_dt_hg38 = readRDS(system.file('extdata', 'repeatmasker_hg38_SINE.bed', package = 'InProgress'))
SINE_dt_hg38 = readRDS("/xchip/beroukhimlab/wolu/testing_svaba/extdata/repeatmasker_hg38_SINE.bed")

#' Line repeatmaskers
#'
#' repeat masker line elements in hg19
#' @name LINE_dt_hg19
#' @docType data
#' @keywords data internal
#' @format \code{data.table}
#LINE_dt_hg19 = readRDS(system.file('extdata', 'repeat_masker_hg19_LINE.bed', package = 'InProgress'))
LINE_dt_hg19 = readRDS("/xchip/beroukhimlab/wolu/testing_svaba/extdata/repeat_masker_hg19_LINE.bed")

#' Sine repeatmaskers
#'
#' repeat masker sine elements in hg19
#' @name SINE_dt_hg19
#' @docType data
#' @keywords data internal
#' @format \code{data.table}
#SINE_dt_hg19 = readRDS(system.file('extdata', 'repeat_masker_hg19_SINE.bed', package = 'InProgress'))
SINE_dt_hg19 = readRDS("/xchip/beroukhimlab/wolu/testing_svaba/extdata/repeat_masker_hg19_SINE.bed")

#' Check format of bedpe
#' @name check_reformat
#' @title reformat for matching
#' @param i iteration value
#' @param df bedpe to be passed 
#' @return data.table that is ordered with lower position first
#' @description Reorders for accurate distance calculation
#' @import data.table
#' @keywords internal
#' 
check_reformat = function(i, df){
  
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

#' Closest line element
#' @name find_closest_match_line
#' @title annotates the closest line element
#' @param i iteration value
#' @param bedpe_l ordered bedpe from \link[InProgress]{check_reformat}
#' @param LINE_dt line elements
#' @return data.table with closest line element distance
#' @description Annotates the distance to closest line element
#' @import data.table
#' @keywords internal
#' 
find_closest_match_line = function(i, bedpe_l, LINE_dt = NULL, LINE_dt_ranges = NULL){

  row_l <- bedpe_l[i,] #the particular row of the bedpe file we are working with

  
  ### both bedpe and ref should be sorted so lower bkpt comes first 
  ref_sub <- LINE_dt[ seqnames == row_l$chrom1 & seqnames == row_l$chrom2]
  ref_sub_ranges <- LINE_dt_ranges[seqnames(LINE_dt_ranges) == row_l$chrom1 | seqnames(LINE_dt_ranges) == row_l$chrom2]
  #ref_sub_ranges <- makeGRangesFromDataFrame(ref_sub, seqnames.field="seqnames", start.field="start", end.field ="end",strand.field="V4") changed this line to increase efficiency

  row_l_str1 <- GRanges(row_l$chrom1, IRanges(as.numeric(as.character(row_l$start1)),width=1))
  row_l_str2 <- GRanges(row_l$chrom2, IRanges(as.numeric(as.character(row_l$start2)),width=1))
  
  check_str1_overlap <- ref_sub_ranges %&% row_l_str1
  check_str2_overlap <- ref_sub_ranges %&% row_l_str2
  
  ref_sub[,str1_dist_l := ifelse(length(check_str1_overlap)>=1, 0, min(abs(start - as.numeric(as.character(row_l$start1))), abs(end - as.numeric(as.character(row_l$start1)))))]
  ref_sub[,str2_dist_l := ifelse(length(check_str2_overlap)>=1, 0, min(abs(start - as.numeric(as.character(row_l$start2))), abs(end - as.numeric(as.character(row_l$start2)))))]
  
  str1_min <- min(ref_sub$str1_dist_l)
  str2_min <- min(ref_sub$str2_dist_l)
  line_dist <- str1_min+str2_min
  if(line_dist == Inf){ #modeling the case where there is no minimum;both str1_dist and str2_dist columns will be empty
    return(cbind(row_l, line_dist=""))
  } else {
    return (cbind(row_l, line_dist))
  }
}

#' Closest sine element
#' @name find_closest_match_sine
#' @title annotates the closest sine element
#' @param i iteration value
#' @param bedpe_s ordered bedpe from \link[InProgress]{check_reformat}
#' @param SINE_dt sine elements
#' @return data.table with closest sine element distance
#' @description Annotates the distance to closest sine element
#' @import data.table
#' @keywords internal
#' 
find_closest_match_sine = function(i, bedpe_s, SINE_dt=NULL, SINE_dt_ranges = NULL){
  row_s <- bedpe_s[i,] 
  
  ### both bedpe and ref should be sorted so lower bkpt comes first 
  ref_sub <- SINE_dt[ seqnames == row_s$chrom1 & seqnames == row_s$chrom2]
  ref_sub_ranges <- SINE_dt_ranges[seqnames(SINE_dt_ranges) == row_s$chrom1 | seqnames(SINE_dt_ranges) == row_s$chrom2]
  #ref_sub_ranges <- makeGRangesFromDataFrame(ref_sub, seqnames.field="seqnames", start.field="start", end.field ="end",strand.field="V4")
  
  row_s_str1 <- GRanges(row_s$chrom1, IRanges(as.numeric(as.character(row_s$start1)),width=1))
  row_s_str2 <- GRanges(row_s$chrom2, IRanges(as.numeric(as.character(row_s$start2)),width=1))
  
  check_str1_overlap <- ref_sub_ranges %&% row_s_str1
  check_str2_overlap <- ref_sub_ranges %&% row_s_str2
  
  ref_sub[,str1_dist_s := ifelse(length(check_str1_overlap)>=1, 0, min(abs(start - as.numeric(as.character(row_s$start1))), abs(end - as.numeric(as.character(row_s$start1)))))]
  ref_sub[,str2_dist_s := ifelse(length(check_str2_overlap)>=1, 0, min(abs(start - as.numeric(as.character(row_s$start2))), abs(end - as.numeric(as.character(row_s$start2)))))]
  
  str1_min <- min(ref_sub$str1_dist_s)
  str2_min <- min(ref_sub$str2_dist_s)
  sine_dist <- str1_min+str2_min
  if(sine_dist == Inf){ #modeling the case where there is no minimum;both str1_dist and str2_dist columns will be empty
    return(cbind(row_s, sine_dist=""))
  }
  else {
    return (cbind(row_s, sine_dist))
  }
}

#' Wrapper for LINE/SINE annotation
#' @name closest_line_sine 
#' @title Annotate LINE and SINE elements
#' @param bp \href{https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format}{Bedpe} from \link[InProgress]{svaba_vcf2bedpe} or \link[InProgress]{manta_vcf2bedpe}
#' @param cores Number of cores to run on, default is 1
#' @param genome run under hg19 or hg38
#' @return \href{https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format}{Bedpe} with a columns added for distance to nearest LINE and SINE element
#' @description
#' 
#' Annotates the distance to nearest LINE and SINE element using repeat masker elements from \href{http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html}{RepeatMasker}.
#' Annotation is in hg38.
#' 
#' @import data.table
#' @importFrom parallel mclapply
#' @export
closest_line_sine = function(bp = NULL, genome = NULL, cores = 1) {
  
  if(is.null(bp)) {
    stop('NULL input')
  }
  
  if(as.character(genome) == 'hg19') {
    LINE_dt = LINE_dt_hg19
    colnames(LINE_dt)[1:3] <- c('seqnames','start','end')
    LINE_dt_ranges = makeGRangesFromDataFrame(LINE_dt, seqnames.field="seqnames", start.field="start", end.field ="end")
    SINE_dt = SINE_dt_hg19
    colnames(SINE_dt)[1:3] <- c('seqnames','start','end')
    SINE_dt_ranges = makeGRangesFromDataFrame(SINE_dt, seqnames.field="seqnames", start.field="start", end.field ="end")
  } else if(as.character(genome) == 'hg38') {
    LINE_dt = LINE_dt_hg38
    LINE_dt_ranges = makeGRangesFromDataFrame(LINE_dt, seqnames.field="seqnames", start.field="start", end.field ="end")
    SINE_dt = SINE_dt_hg38
    SINE_dt_ranges = makeGRangesFromDataFrame(SINE_dt, seqnames.field="seqnames", start.field="start", end.field ="end")
  } else {
    stop('Please state hg19 or hg38 as genome')
  }
  
  cat("Checking format...")
  bp_ord <- rbindlist(mclapply(1:nrow(bp), check_reformat, bp, mc.cores = cores))
  cat("done.\n")
  bp_ord[chrom1 == 23, chrom1 := "X"]
  bp_ord[chrom2 == 23, chrom2 := "X"]
  bp_ord[chrom1 == 24, chrom1 := "Y"]
  bp_ord[chrom2 == 24, chrom2 := "Y"]
  cat("Comparing against LINE elements...")
  line_annotated <- rbindlist(mclapply(1:nrow(bp_ord), find_closest_match_line, bp_ord, LINE_dt = LINE_dt, LINE_dt_ranges = LINE_dt_ranges, mc.cores = cores))
  cat("done.\n")
  
  cat("Comparing against SINE elements...")
  sine_line_annotated <- rbindlist(mclapply(1:nrow(line_annotated), find_closest_match_sine, line_annotated, SINE_dt = SINE_dt, SINE_dt_ranges = SINE_dt_ranges, mc.cores = cores))
  cat("done.\n")
  
  return(sine_line_annotated)
}



run_it <- function(lof_pth) {
  sample.name <- unlist(strsplit(lof_pth, "filtered_only/"))[2]
  print(sample.name) #added this line 
  output_path <- paste0('/xchip/beroukhimlab/wolu/testing_svaba/outputs/final_fuzzymatch/fuzzymatch_only/',sample.name,'_fuzzy.bedpe')
  
  bedpe <- fread(lof_pth)
  bedpe[, chrom1 := gsub('chr','',chrom1)]
  bedpe[, chrom2 := gsub('chr','',chrom2)]
  
  fuzzy <- closest_germline(bp = bedpe, cores = 5, genome = 'hg19')
  line_sine <- closest_line_sine(bp = fuzzy, genome = 'hg19', cores = 5) 
  
  write.table(line_sine, output_path, sep = '\t', row.names = F, col.names = T, quote = F)
}

filepath <- commandArgs(T)[1]
run_it(filepath)