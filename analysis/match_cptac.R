library(data.table)

#cptac_dir<-'/Volumes/xchip_beroukhimlab/Simona/SVsigs/data/'
#manifest<-fread(paste0(cptac_dir, 'gdc_manifest_20220707_181137.txt'))
#metadata<-fread(paste0(cptac_dir, 'gdc_sample_sheet.2022-07-22.tsv'))

manifest<-fread('/Users/shu/germline_svm/data/gdc_manifest_20220707_181137.txt')
metadata<-fread('/Volumes/xchip_beroukhimlab/Simona/SVsigs/data/gdc_sample_sheet.2022-08-02.tsv')



single_case<-sapply(metadata$`Case ID`, function(i){
  split_case<-strsplit(i, ",")
  return(split_case[[1]][1])
})
metadata$single_case<-single_case

#fc-secure-dc1febc1-419d-4083-882c-5a6058c1edef/cptac/data/9a5c7c77-02ed-4bf2-bf16-f08947fabadd/35cbfe6b-e6d9-4191-8f9f-c1a53233e611_gdc_realn.bai

#files are in folder like so

#gs://fc-secure-dc1febc1-419d-4083-882c-5a6058c1edef/cptac/data/id(from manifest)/filename(from manifest)

#gsutil cp gs://fc-secure-dc1febc1-419d-4083-882c-5a6058c1edef/cptac/data/035ca893-49e4-4b59-8f08-ea13e64c7b9c/74bafef4-c6ea-4079-94d1-272f647dc377_targetedsequencing_gdc_realn.bam .
gs_path<-'gs://fc-secure-dc1febc1-419d-4083-882c-5a6058c1edef/cptac/data/'

#NOTE TARGETED SEQUENCING FILES ARE CURRENTLY NOT INDEXED ONLY BAM NO BAI 

match_cptac<-function(case_id){
  case_df<-metadata[metadata$single_case==case_id]
  test<-as.character(case_df$`Sample Type`)
  #print(test)
  if ( any(test %like% "Tumor") & any(test %like% "Blood Derived Normal")){
    
    tumor_idx<-which(test %like% "Tumor")
    normal_idx<-which(test %like% "Blood Derived Normal")

    tumor_df<-case_df[tumor_idx,]
    #print(tumor_df)
    panel_row<-tumor_df[which(tumor_df$`File Name` %like% "targetedsequencing"),]
    
    tumor_wgs_idx<-which(tumor_df$`File Name` %like% "wgs_gdc_realn")
    #if multiple tumors from sample, just pick one 
    #print(tumor_wgs_idx)
    #print(tumor_wgs_idx[1])
    
    tumor_row<-tumor_df[tumor_wgs_idx,]

    normal_row<-case_df[normal_idx,]
    
    print(tumor_row$`File Name`)
    
    tumor_manifest_id<-manifest$id[manifest$filename==tumor_row$`File Name`]
    panel_manifest_id<-manifest$id[manifest$filename==panel_row$`File Name`]
    normal_manifest_id<-manifest$id[manifest$filename==normal_row$`File Name`]
    
    
    normal_file_id<-manifest$filename[manifest$filename==tumor_row$`File Name`]
    tumor_file_id<-manifest$filename[manifest$filename==tumor_row$`File Name`]
    panel_file_id<-manifest$filename[manifest$filename==tumor_row$`File Name`]
    

    return(data.table(case_id, 
          normal_sample_id=normal_row$`Sample ID`, normal_manifest_id, normal_file_id, 
          tumor_sample_id=tumor_row$`Sample ID`, tumor_manifest_id, tumor_file_id, 
          panel_sample_id=panel_row$`Sample ID`, panel_manifest_id, panel_file_id, fill=T))
  }
  
  else{
    print (paste0(case_id, " has no tumor normal pair"))
  }
  
}


cptac_pair<-lapply(unique(metadata$single_case), match_cptac)
cptac_pair_rbind<-rbindlist(cptac_pair)


#THERE ARE SOME DUPLICATES
#some samples have multiple wgs primary tumor samples, just pick one
cptac_pair_df<-cptac_pair_rbind[!duplicated(cptac_pair_rbind[,c('case_id','normal_sample_id')]),]





#path for normals
cptac_pair_df[, PATH_normal_bam:=paste0(gs_path, normal_manifest_id, '/', normal_file_id)]
cptac_pair_df[, PATH_normal_bai:=paste0(substr(PATH_normal_bam, 1, nchar(PATH_normal_bam)-4), '.bai')]
#path for wgs tumors
cptac_pair_df[, PATH_tumor_bam:=paste0(gs_path, tumor_manifest_id, '/', tumor_file_id)]
cptac_pair_df[, PATH_tumor_bai:=paste0(substr(PATH_tumor_bam, 1, nchar(PATH_tumor_bam)-4), '.bai')]
#panel doesn't have bai, needs to be indexed 
cptac_pair_df[, PATH_panel_bam:=paste0(gs_path, panel_manifest_id, '/', panel_file_id)]

setnames(cptac_pair_df, "case_id", "entity:sample_id")
#write to tsv to upload to terra 

#fwrite(cptac_pair_df, '/Users/shu/germline_svm/data/20220801_cptac_table.tsv', sep='\t')


#get subset of cptac that has both panel and tumor

cptac_pair_df_subset<-cptac_pair_df[(!is.na(cptac_pair_df$tumor_sample_id)) & (!is.na(cptac_pair_df$panel_sample_id))]
#fwrite(cptac_pair_df_subset, file='/Users/shu/germline_svm/data/20220801_cptac_table_subset_panelandtumor.tsv', sep='\t')

#


#read in complete file for now
test<-fread('/Users/shu/germline_svm/data/20220801_cptac_table.tsv')


