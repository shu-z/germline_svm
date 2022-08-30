library(data.table)

#########################################
# code to filter events >1k (don't need to run every time)
# aims to select equal # of SVs from equal # of tumor types 


#read in pcawg metadata file
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

#this is updated one as of 04/19/2022, with updated coding annotations
#has extra filtering, and cn/exon annotations
snowman_germdist_dir<-'/Volumes/xchip_beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/data/4_annot/'
all_paths<-paste0(snowman_germdist_dir, list.files(snowman_germdist_dir))
map_df<-rbindlist(lapply(all_paths, map_id), fill=T)
sort(table(map_df$project_code) )  #min num tumors is 7 from dlbcl

#write intermediate mapping file so don't need to read again 
#write.csv(map_df, '/Volumes/xchip_beroukhimlab/Shu/ccle/20220425_snowman_mapping_newannot.csv')
map_df<-fread('/Volumes/xchip_beroukhimlab/Shu/ccle/20220425_snowman_mapping_newannot.csv', header=T, drop=1)

#select sample number of files per tumor type 
even_tt_sampling<-function(row, test_train){
  tumor_type<-row['V1']
  ifelse(test_train=='test', n_samp<-row['n_test'], n_samp<-row['n_train'])
  tumor_paths<-map_df$path[map_df$project_code == tumor_type]
  return(data.table(paths=sample(tumor_paths, size=n_samp), project_code=tumor_type))
}

#select tumors for test set
n_samp_tt<-data.table(sort(table(map_df$project_code)))  #min num tumors is 7 from dlbcl 
n_samp_tt[, n_test:=ifelse(N<15, N, ifelse(N<25, N-15, 10))]  
n_samp_tt[, n_train:=ifelse(N<15, 0, 15)]  


#this will be consistent test set for following analyses as of 20220601
test_samples<-rbindlist(apply(n_samp_tt, 1, even_tt_sampling, 'test'))
#write.csv(test_samples, '/Volumes/xchip_beroukhimlab/Shu/svm/data/20220601_snowman_newannot_test_tumorslist.csv')

#remove test samples from map_df
map_df_train<-map_df[!(map_df$path %in% test_samples$paths)]

#randomly select tumors for train set (all train SVs will sample from these tumors)
train_samples<-rbindlist(apply(n_samp_tt, 1, even_tt_sampling, test_train='train'))
train_samples<-train_samples[!(is.na(train_samples$paths))]


#function for selecting SVs
#remove short events and select same number of svs per tumor
even_sv_sampling<-function(row, n_max){
  sample_df<-fread(row[1])
  long_sample_df<-sample_df[(sample_df$SPAN>=1e3) | (sample_df$SPAN ==-1)]
  print(nrow(long_sample_df))
  n_len<-min(nrow(long_sample_df), n_max)
  return(cbind(long_sample_df[sample(seq_len(nrow(long_sample_df)), size = n_len),], 
               n_long_events=nrow(long_sample_df), project_code=row[2], file_path=row[1]))
}


#get all SVs from test tumors, but sample same # from train tumors 
test_samples_df<-rbindlist(apply(test_samples, 1, even_sv_sampling, n_max=100000))
train_samples_df<-rbindlist(apply(train_samples, 1, even_sv_sampling, n_max=600))

#write.csv(test_samples_df, '/Volumes/xchip_beroukhimlab/Shu/svm/data/20220601_snowman_newannot_test_10tumors_allsvs.csv')
#write.csv(train_samples_df, '/Volumes/xchip_beroukhimlab/Shu/svm/data/20220601_snowman_newannot_train_15tumors_600svs.csv')






################### notes to self

#throw out dlbcl
#sample from same number still
#test on dlbcl
#find new median of SVs to sample (613)
#look into tumor type representation in somatic total df
#look at tumor type aucs with differnet num train svs vs global auc
# same test set sv
#run multiple times
