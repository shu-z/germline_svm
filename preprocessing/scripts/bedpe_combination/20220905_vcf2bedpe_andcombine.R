require(data.table)
require(parallel)

##### functioons #####
vcf_to_dt <- function(vcf_path) {
  cat(paste0(vcf_path, "\n"))
  if (!file.exists(vcf_path)) {
    print(paste("File does not exist",vcf_path))
  }
  
  # Don't continue if empty GRanges object...
  #if (!as.numeric(system(paste("grep -v '^#'", vcf_path, "| wc -l"), intern=TRUE))) {
  
  #  return (GRangesList())
  #}
  
  # Read in vcf as data.table object (vcf_dt)...
  # This grep command returns all lines of the file that don't start with "#"...
  cat("Reading file...\n")
  vcf_dt <- fread(cmd=paste("grep -v '^#'", vcf_path),sep='\t')
  
  # Set colnames of vcf_dt to standard...
  if (nrow(vcf_dt) == 0) {
    return (vcf_dt)
  }
  if (ncol(vcf_dt)==10) {
    setnames(vcf_dt, paste0("V",seq(1:10)), c("seqnames","start","ID","REF","ALT","QUAL","FILTER","INFO","GENO","NORMAL"))
  } else {
    setnames(vcf_dt, paste0("V",seq(1:11)), c("seqnames","start","ID","REF","ALT","QUAL","FILTER","INFO","GENO","NORMAL","TUMOR"), skip_absent=TRUE)
  }
  
  cat("Gathering Metadata...\n")
  if ("INFO" %in% colnames(vcf_dt) ) {
    # Extract info from columns... 
    vcf_dt[, SPAN := as.numeric(gsub(".*?SPAN=([-0-9]+).*","\\1",INFO))]
    vcf_dt$sample = gsub("(.*?)_.*","\\1",basename(vcf_path))
    vcf_dt[, uid := gsub("([0-9]+):(1|2)", "\\1", ID)]
    vcf_dt[, EVDNC := gsub(".*?EVDNC=([A-Z]+).*", "\\1", INFO)]
    vcf_dt[, NUMPARTS := as.integer(gsub(".*?NUMPARTS=([0-9]+).*", "\\1", INFO))]
    vcf_dt[, SCTG := gsub(".*?SCTG=(.*?);.*", "\\1", INFO)]
    vcf_dt[, MAPQ := as.integer(gsub(".*?;MAPQ=([0-9]+).*", "\\1", INFO))]
    vcf_dt[, REPSEQ := gsub(".*?;REPSEQ=([A-Z]+).*", "\\1", INFO)]
    vcf_dt[, REPSEQ := ifelse(grepl(";", REPSEQ), "", REPSEQ)] 
    vcf_dt[, HOMSEQ := gsub(".*?;HOMSEQ=([A-Z]+).*", "\\1", INFO)] 
    vcf_dt[, HOMSEQ := ifelse(grepl(";", HOMSEQ), "", HOMSEQ)] 
    vcf_dt[, INSERTION := gsub(".*?;INSERTION=([A-Z]+).*", "\\1", INFO)] 
    vcf_dt[, INSERTION := ifelse(grepl(";", INSERTION), "", INSERTION)]
    vcf_dt[, NDISC := as.numeric(gsub(".*?NDISC=([0-9]+).*", "\\1", INFO))]
    vcf_dt[, SVMETHOD := substr(INFO,regexpr("SVMETHOD=",INFO)+nchar("SVMETHOD="),regexpr(";NDISC",INFO)-1)]
    
  }
  
  # More extraction regexpr stuff...
  if ("TUMOR" %in% colnames(vcf_dt)) {
    vcf_dt[, TUMALT :=  as.integer(strsplit(TUMOR, ":")[[1]][2]) , by=uid]
    vcf_dt[, TUMCOV :=  as.integer(strsplit(TUMOR, ":")[[1]][3]) , by=uid]
    vcf_dt[, TUMLOD :=  as.numeric(strsplit(TUMOR, ":")[[1]][9]) , by=uid]
  }
  if ("NORMAL" %in% colnames(vcf_dt)) {
    vcf_dt[, NORMCOV :=  as.integer(strsplit(NORMAL, ":")[[1]][3]) , by=uid]
    vcf_dt[, NORMALT :=  as.integer(strsplit(NORMAL, ":")[[1]][2]) , by=uid]
    vcf_dt[, NORMLOD :=  as.numeric(strsplit(NORMAL, ":")[[1]][9]) , by=uid]
  }
  
  cat("Cleaning up...\n")
  vcf_dt[, strand := ifelse(grepl("^\\[", ALT) | grepl("^\\]", ALT), '-', '+')]
  vcf_dt[, inv := strand[1] == strand[2], by=uid]
  vcf_dt[, altstrand := rev(strand), by=uid]
  vcf_dt[, altpos := as.integer(gsub(".*?:([0-9]+).*", "\\1", ALT))]
  vcf_dt[, altchr := gsub(".*?(\\[|\\])(.*?):([0-9]+).*", "\\2", ALT)]
  vcf_dt[, end := start] #why are we making the start and end locations of the chromosomes the same?
  
  bad.ix <- vcf_dt[grepl("^G|^M", seqnames), uid]
  vcf_dt <- vcf_dt[!uid %in% bad.ix]
  
  # Set SID == name of file path without the extension...
  vcf_dt[, sid := basename(tools::file_path_sans_ext(vcf_path))]
  
  # add chr prefix
  vcf_dt$seqnames <- paste0("chr",vcf_dt$seqnames)
  #adding a print statement here to see if the code goes beyond the rough input
  #print(vcf_dt)
  return(vcf_dt)
}

build_bedpe_with_metadata <- function(merged_dt) {
  cat("Building bedpe...\n")
  
  ### get mate indexes
  merged_dt[, mates_idx := unlist(strsplit(ID, ":"))[1], by = "uid"]
  #changing the argument that goes into which_mate from "3" too "2" because there is only one : in the ID
  #print(merged_dt)
  #iterlast copying this to the for loop
  merged_dt[, which_mate := unlist(strsplit(ID, ":"))[2], by = "uid"]
  #print(c(merged_dt$which_mate),c(merged_dt$seqnames)) #added this print statement to see what this is doing
  temp_bedpe <- NULL
  removed_bnd <- NULL
  for(i in 1:length(unique(merged_dt$mates_idx))){
    
    foo <- merged_dt[mates_idx == unique(merged_dt$mates_idx)[i]]
    #print(foo) #added this print statement to see what foo looks like before the modifications below
    #print(nrow(foo)) #adding this print to see what the number of rows that is checked in the if statement is
    
    if(!(nrow(foo)== 2)) {
      mes <- paste0("Breakpoint ",  unique(merged_dt$mates_idx)[i], " has incorrect number of mates for ", foo$sample, " It has been removed.")
      # excludes this breakpoint from 
      continue = FALSE
      removed_bnd <- rbind(removed_bnd, foo)
      warning(mes[1])
    }
    else {
      continue = TRUE
    }
    if(continue) {
      #### build bedpe
      #moved here to update the way foo is derived
      foo[,split_ID := c(1:length(foo$seqnames))] #adding this so that there are unique values with which which_mate is derived; adds a lot of complexity though
      foo[, which_mate := unlist(strsplit(ID, ":"))[2], by = "split_ID"] #doing it twice so that it maintains the 2 paired mates but updates the which_mate variable
      #print(foo) #want to see how foo looks with the secondary modification for the which_mates variable;unfortunately, if a chr has a mate that is itself, then it defaults to
      #mark as the same mate. We need some other property that is unique to the mates to tag them
      foo1 <- foo[which_mate == 1]
      #print(foo1) #added this print statement
      foo2 <- foo[which_mate == 2]
      #print(foo2) #added this print statement
      
      bedpe_base <- as.data.frame(cbind(foo1$seqnames, foo1$start, foo1$end,
                                        foo2$seqnames, foo2$start, foo2$end))
      #print(bedpe_base) #added this print statement 'cause the code still insists that one of these is empty
      colnames(bedpe_base) <- c("chrom1", "start1", "end1", "chrom2","start2","end2")
      
      if(!(foo1$sid == foo2$sid)){
        stop("Multiple samples are being processed, one at a time please...")
      }
      
      bedpe_base <- cbind(bedpe_base, paste0(foo$sample[1],"_", foo$mates_idx[1]))
      colnames(bedpe_base)[7] <- "name"
      
      bedpe_base <- cbind(bedpe_base, foo$QUAL[1])
      colnames(bedpe_base)[8] <- "score"
      
      bedpe_base <- cbind(bedpe_base, foo1$strand[1])
      bedpe_base <- cbind(bedpe_base, foo2$strand[1])
      colnames(bedpe_base)[9:10] <- c("strand1", "strand2")
      
      refs_alts <- as.data.frame(cbind(foo1$REF[1], foo1$ALT[1],
                                       foo2$REF[1], foo2$ALT[1]))
      colnames(refs_alts) <- c("REF_1","ALT_1","REF_2","ALT_2")
      bedpe_base <- cbind(bedpe_base, refs_alts)
      bedpe_base <- cbind(bedpe_base, foo1[,c("SPAN", "HOMSEQ","INSERTION","NDISC","FILTER","sample", "TUMALT", 'INFO', 'GENO','TUMOR')])
      mapqs <- as.data.frame(cbind(foo1$MAPQ[1], foo2$MAPQ[1]))
      colnames(mapqs) <- c("MAPQ_1","MAPQ_2")
      bedpe_base <- cbind(bedpe_base, mapqs)
      
      
      temp_bedpe <- rbind(temp_bedpe, bedpe_base)
    }
  }
  
  bedpe <- as.data.table(temp_bedpe)
  return(bedpe)
}


merge <- function(lof_pth) {
  path_germline <- unlist(strsplit(lof_pth, "%"))[1]
  path_somatic <- unlist(strsplit(lof_pth, "%"))[2]
  
  somatic_vcf_dt <- vcf_to_dt(path_somatic)
  germline_vcf_dt <- vcf_to_dt(path_germline)
  
  somatic_bedpe <- build_bedpe_with_metadata(somatic_vcf_dt)
  germline_bedepe <- build_bedpe_with_metadata(germline_vcf_dt)
  
  somatic_bedpe[, CLASS := "SOMATIC"]
  germline_bedepe[, CLASS := "GERMLINE"]
  
  test_bed <- rbind(somatic_bedpe, germline_bedepe)
  somatic_id <- unlist(strsplit(unlist(strsplit(path_somatic, "/"))[6], "[.]"))[1] #changing the first subsetting from 9 to 6 because that is where the 'unique' value for each sample lies
  germline_id <-  unlist(strsplit(unlist(strsplit(path_germline, "/"))[6], "[.]"))[1]
  
  # Change to own output directory
  if(somatic_id == germline_id) {
    output <- paste0("/xchip/beroukhimlab/siyun/germline_classifier/outputs/bedpe_combination/", somatic_id, "_combined_germ_soma_test.bedpe")
    
    write.table(test_bed, output, row.names = F, col.names = T, sep = "\t", quote = F)
  }
}


# Takes in the text file of all sample paths
# DONT NEED THIS
merge_multiple <- function(samples_file) {
  
  all_samples <- file(samples_file, open="r")
  lines <- readLines(all_samples)
  for (i in 1:length(lines)){
  
    path_germline <- unlist(strsplit(lines[i], "%"))[1]
    path_somatic <- unlist(strsplit(lines[i], "%"))[2]
    
    somatic_vcf_dt <- vcf_to_dt(path_somatic)
    germline_vcf_dt <- vcf_to_dt(path_germline)
    
    somatic_bedpe <- build_bedpe_with_metadata(somatic_vcf_dt)
    germline_bedepe <- build_bedpe_with_metadata(germline_vcf_dt)
    
    somatic_bedpe[, CLASS := "SOMATIC"]
    germline_bedepe[, CLASS := "GERMLINE"]
    
    if (i == 1) {
      test_bed <- rbind(somatic_bedpe, germline_bedepe)
    }
    else {
      temp_bed <- rbind(somatic_bedpe, germline_bedepe)
      test_bed <- rbind(test_bed, temp_bed)
    }
    
  }
  output <- paste0("/xchip/beroukhimlab/siyun/germline_classifier/outputs/bedpe_combination/all_samples_combined_germ_soma_test.bedpe")
  write.table(test_bed, output, row.names = F, col.names = T, sep = "\t", quote = F)
  
  close(all_samples)
}


# This is the MAIN!
if(length(commandArgs(T) > 0)) {
  file = commandArgs(T)[1]
  all_samples <- file(file, open="r")
  lines <- readLines(all_samples)
  
  for (i in 1:length(lines)){
    merge(lines[i])
  }
  
} else {
  stop("Need file path")
}

