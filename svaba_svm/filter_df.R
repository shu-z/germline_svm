library(data.table)

#########################################
# code to filter events >1k (don't need to run every time)
# aims to select equal # of SVs from equal # of tumor types 


#read in pcawg metadata file
svaba_metadata <- fread("Z:/siyun/testing/pcawg_svsigs_metadata.csv",sep=',', header = TRUE) #germline_metadata because it should be a lot more comprehensive
colnames(svaba_metadata) <- c('germline_path','somatic_path','Abbreviated','TCGA_ID','Tumor_type')



#function to map sample name to type 
map_id<-function(i){
  sample_df<-fread(i, nrows=1)
  sample_vcf_path <- paste0("/xchip/beroukhimlab/Shu/tcga_svaba_germline_marcin/",unlist(strsplit(sample_df$sample,"[.]"))[1],".svaba.germline.sv.vcf")
  sample_name <- paste0(unlist(strsplit(sample_df$sample,"[.]"))[1],".svaba.germline.sv.vcf")
  tumor_type <- svaba_metadata$Tumor_type[svaba_metadata$germline_path %in% sample_vcf_path]
  return(data.table(cbind(path=i, sample_name, sample_vcf_path, tumor_type)))
}

#this is updated one as of 04/19/2022, with updated coding annotations;comment from snowman_filter
#has extra filtering, and cn/exon annotations
svaba_annotations_dir <- "Z:/wolu/testing_svaba/outputs/final_annot/annotated_only/"
all_paths<-paste0(svaba_annotations_dir, list.files(svaba_annotations_dir))
map_df<-rbindlist(lapply(all_paths, map_id), fill=T)
sort(table(map_df$tumor_type))  #min num tumors is 7 from dlbcl

#write intermediate mapping file so don't need to read again 
#write.csv(map_df,"Z:/wolu/testing_svaba/outputs/svm/20221026_svaba_mapping_newannot.csv")

map_df <- fread("Z:/wolu/testing_svaba/outputs/svm/20221026_svaba_mapping_newannot.csv",header = TRUE, drop = 1)

#select sample number of files per tumor type 
even_tt_sampling<-function(row, test_train){ #this function might be useless but let's figure out why
  tumor_type<-row['V1']
  ifelse(test_train=='test', n_samp<-row['n_test'], n_samp<-row['n_train'])
  tumor_paths<-map_df$path[map_df$tumor_type == tumor_type]
  return(data.table(paths=sample(tumor_paths, size=n_samp), project_code=tumor_type))
}

#select tumors for test set
n_samp_tt<-data.table(sort(table(map_df$tumor_type)))  #min num tumors is 7 from dlbcl
n_samp_tt[, n_test:=ifelse(N<15, N, ifelse(N<25, N-15, 10))]  
n_samp_tt[, n_train:=ifelse(N<15, 0, 15)]  


#this will be consistent test set for following analyses as of 20220601 (comment from previous snowman svm)
test_samples<-rbindlist(apply(n_samp_tt, 1, even_tt_sampling, 'test'))
write.csv(test_samples,"Z:/wolu/testing_svaba/outputs/svm/20221027_svaba_annot_test_tumorlist.csv")

#remove test samples from map_df
map_df_train<-map_df[!(map_df$path %in% test_samples$paths)]

#randomly select tumors for train set (all train SVs will sample from these tumors)
train_samples<-rbindlist(apply(n_samp_tt, 1, even_tt_sampling, test_train='train'))
train_samples<-train_samples[!(is.na(train_samples$paths))]
write.csv(train_samples, "Z:/wolu/testing_svaba/outputs/svm/20221031_svaba_annot_train_tumorlist.csv")


#function for selecting SVs
#remove short events and select same number of svs per tumor
even_sv_sampling<-function(row, n_max){
  sample_df<-fread(row[1])
  long_sample_df<-sample_df[(sample_df$SPAN>=1e3) | (sample_df$SPAN ==-1)] #select the SVs for which the span is >1000bp or a translocation from a sample
  print(nrow(long_sample_df)) #how many Svs are we getting from that sample
  n_len<-min(nrow(long_sample_df), n_max)
  return(cbind(long_sample_df[sample(seq_len(nrow(long_sample_df)), size = n_len),], 
               n_long_events=nrow(long_sample_df), project_code=row[2], file_path=row[1]))
}


#get all SVs from test tumors, but sample same # from train tumors 
test_samples_df<-rbindlist(apply(test_samples, 1, even_sv_sampling, n_max=100000))#this assumes that a single file doesn't have up to 100000 SVs;
#this is true because the max no. of SV between germ and somatic events is 35795 (order of 10^4)
train_samples_df<-rbindlist(apply(train_samples, 1, even_sv_sampling, n_max=700))
#the cutoff is selected based on the median no. of SVs amongst the train set to mitigate overrepresentation



write.csv(test_samples_df, 'Z:/wolu/testing_svaba/outputs/svm/20221027_svaba_annot_test_10tumors_allsvs.csv')
write.csv(train_samples_df, 'Z:/wolu/testing_svaba/outputs/svm/20221027_svaba_annot_train_15tumors_700svs.csv')






################### notes from Shu's snowman implementation.

#throw out dlbcl
#sample from same number still
#test on dlbcl
#find new median of SVs to sample (613)
#look into tumor type representation in somatic total df
#look at tumor type aucs with differnet num train svs vs global auc
# same test set sv
#run multiple times