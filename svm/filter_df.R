#########################################
# code to filter events >1k (don't need to run every time)
# aims to select equal # of SVs from equal # of tumor types 

pcawg_metadata<-fread('/Volumes/xchip_beroukhimlab/Alex/PCAWG/PCAWG_metadata_OCT2020.csv', sep=',')
#pcawg_metadata<-fread('/Users/shu/ClusterSV/supp_table_pcawg.csv', sep=',')
#which(pcawg_metadata == '00493087-9d9d-40ca-86d5-936f1b951c93', arr.ind=T)

#function to map sample name to type 
map_id<-function(i){
  sample_df<-fread(i, nrows=1)
  sample_name<-substring(sample_df$sample, 1, 36)
  project_code<-pcawg_metadata$dcc_project_code[pcawg_metadata$tumor_wgs_submitter_sample_id %in% sample_name]
  return(data.table(cbind(path=i, sample_name, project_code)))
}

#read in snowman germline somatic calls 
#new one has line-sine annotation
snowman_germdist_dir<-'/Volumes/xchip_beroukhimlab/Alex/CCLE_PCAWG/Germline_Correctiion/Test_Set/Germ_dist/FUNCT_ANNOT/Line_Sine_Annot/'
all_paths<-paste0(snowman_germdist_dir, list.files(snowman_germdist_dir))
map_df<-rbindlist(lapply(all_paths, map_id), fil=T)
#table(map_df$project_code)   #min num tumors is 7 from dlbcl 

#select sample number of files per tutmor type 
even_tt_sampling<-function(tumor_type){
  tumor_paths<-map_df$path[map_df$project_code == tumor_type]
  return(data.table(paths=sample(tumor_paths, size=7), project_code=tumor_type))
}

#remove short events and select same number of svs 
even_sv_sampling<-function(row){
  sample_df<-fread(row[1])
  long_sample_df<-sample_df[(sample_df$SPAN>=1e3) | (sample_df$SPAN ==-1)]
  print(nrow(long_sample_df))
  n_len<-min(nrow(long_sample_df), 1200)
  return(cbind(long_sample_df[sample(seq_len(nrow(long_sample_df)), size = n_len),], 
               n_long_events=nrow(long_sample_df), project_code=row[2], file_path=row[1]))
}

sample_paths<-rbindlist(lapply(unique(map_df$project_code), even_tt_sampling))
final_svm_df<-rbindlist(apply(sample_paths, 1, even_sv_sampling))
write.csv(final_svm_df, '/Volumes/xchip_beroukhimlab/Shu/ccle/20220302_snowman_evenlysampledtumor.csv')

