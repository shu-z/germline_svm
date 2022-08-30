library(data.table)
library(e1071)
library(caTools)
library(ROCR)
library(stringr)
library(ggfortify)
library(caret)
library(ggplot2)


#####script to see if training and testing on different ancestry groups works 

#gnomad SVs do not have ancestry data 
#gnomad_sv<-fread('/Users/shu/Downloads/gnomad_v2.1_sv.sites.bed', nrows=10)


#read in metadata and ancestry data 
pcawg_metadata<-fread('/Volumes/xchip_beroukhimlab/Alex/PCAWG/PCAWG_metadata_OCT2020.csv', sep=',')
ancestry_data<-fread('/Volumes/xchip_beroukhimlab/Alex/CCLE_PCAWG/Germline_Correctiion/Ancestry/pcawg_ancestry_proportions.txt')
ancestry_us<-ancestry_data[ancestry_data$Country=='US']
rm(ancestry_data)


#map ancestry ids to tumor names 
map_ancestry<-apply(ancestry_us, 1, function(row){
  vcf_header<-row['VCF header']
  tumor_wgs_id<-pcawg_metadata$tumor_wgs_submitter_sample_id[pcawg_metadata$normal_wgs_submitter_sample_id %in% vcf_header]
  return(tumor_wgs_id)
})
ancestry_us[, tumor_wgs_id:=map_ancestry]



###############
#use after add_features, once somatic and germline subsets are made 

#df_train_all is all somatic SVs from train set, with equal number of germline SVs
#with 20220601 sample, 62,000 SVs total 
somatic_train<-somatic_subset[somatic_subset$train_test=='train']
germline_train<-germline_subset[germline_subset$train_test=='train']
df_train_all<-rbind(somatic_train, germline_train[sample(seq_len(nrow(germline_train)), size = nrow(somatic_train)),])



features_toscale<-c('homlen', 'insertion_len', 'SPAN', 'gnomad_dist', 'del', 'dup', 'inv', 'inter',
                    'line_dist', 'sine_dist', 'num_sv_sample', 'CN_annot', 'exon_annot')

somatic_test<-somatic_subset[somatic_subset$train_test=='test']
germline_test<-germline_subset[germline_subset$train_test=='test']
df_test_all<-rbind(somatic_test, germline_test)



###########################
test_ancestry<-function(ancestry, classifier){
  ancestry_idx<-(which(ancestry_us[, ..ancestry] >= ancest_cutoff))
  ancestry_subset<-(ancestry_us[ancestry_idx])
  
  test_ancest<-df_test_all[df_test_all$sample_name %in% ancestry_subset$tumor_wgs_id]
  
  test_scaled<-(test_ancest[, lapply(.SD, scale), .SDcols = features_toscale])
  test_scaled<-cbind(test_scaled, sv_class=test_ancest$sv_class)
  
  
  #predict test set with classifier
  y_pred_all <- predict(classifier, newdata = test_scaled, decision.values = T, probability = T)
  probabilities_all<-data.table(attr(y_pred_all, 'probabilities'))
  setcolorder(probabilities_all, c('0', '1'))
  pr_all<-prediction(as.numeric(probabilities_all$`1`), as.numeric(test_scaled$sv_class))
  auc_all <- performance(pr_all, measure = "auc")
  return(data.table(test_ancestry=ancestry, n_test=nrow(test_scaled), auc_all=auc_all@y.values[[1]]))
  
  
}

train_diff_ancestry<-function(ancestry){
  print(paste0('training ancestry: ', ancestry))
  
  ancestry_idx<-(which(ancestry_us[, ..ancestry] >= ancest_cutoff))
  ancestry_subset<-(ancestry_us[ancestry_idx])
  
  df_train_ances<-df_train_all[df_train_all$sample_name %in% ancestry_subset$tumor_wgs_id]
  train <- df_train_ances[sample(seq_len(nrow(df_train_ances)), size = min(n_train, nrow(df_train_ances))), ]

  #scale features in train and test set
  train_scaled<-(train[, lapply(.SD, scale), .SDcols = features_toscale])
  train_scaled<-cbind(train_scaled, sv_class=train$sv_class)


  classifier = svm(formula = sv_class ~ .,
                   data = train_scaled, type = 'C-classification',
                   kernel = 'linear', cachesize=400,
                   probability=T)



  test_results<-rbindlist(lapply(ancestry_all_list, test_ancestry, classifier))
  
  return(cbind(train_ancestry=ancestry, n_train=nrow(train), test_results))
  # test_tt_res<-cbind(test_tt, y_pred_all)
  # split_name<-strsplit(test_tt_res$sample, '[.]')
  # short_name<-lapply(1:length(split_name), function(i){
  #   return(split_name[[i]][1])
  # })
  # test_tt_res$patient_id<-short_name
  # test_tt_res[,missclassified:=ifelse(sv_class==y_pred_all, 0, 1)]

  #get proportion of misclassified SVs per sample s
  #missclass_tt<-rbindlist(lapply(unique(test_tt_res$patient_id), missclass_by_sample, test_df=test_tt_res))




}

#AMR and SAN have 18 and 6 tumors respectively, can test on these but not train 
ancestry_all_list<-c('AFR', 'ASN', 'AMR', 'SAN', 'EUR')
ancestry_train_list<-c('AFR', 'ASN', 'EUR')

n_train<-1000
ancest_cutoff<-0.6
ancestry_svm_df<-rbindlist(lapply(ancestry_all_list, train_diff_ancestry))



#tt_df<-rbindlist(lapply(project_codes[3:5], train_tumor_out))

#write.csv(tt_df, '/Users/shu/germline_svm/data/20220622_train_minusonetumor_testscaled_10ktrain.csv')
#tt_df<-fread('/Users/shu/germline_svm/data/20220622_train_minusonetumor_testscaled_10ktrain.csv')
# tt_unique_auc<-tt_df[!duplicated(tt_df[,c(1,3)]),]
# tt_unique_auc<-tt_unique_auc[,c(1,3)]


#get summary statistics 
tt_med<-aggregate(prop_class ~ tumor_type, tt_df, 'median')
tt_mean<-aggregate(prop_class ~ tumor_type, tt_df, 'mean')
tt_sd<-aggregate(prop_class ~ tumor_type, tt_df, 'sd')
tt_auc<-aggregate(auc_all ~ tumor_type, tt_df, 'median')

tt_summary<-cbind(tt_med, tt_mean[,2], tt_sd[,2], tt_auc[,2])
colnames(tt_summary)<-c('tumor_type', 'median', 'mean', 'sd', 'auc')

######################