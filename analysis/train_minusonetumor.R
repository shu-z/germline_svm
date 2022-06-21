library(data.table)
library(e1071)
library(caTools)
library(ROCR)
library(stringr)
library(ggfortify)
library(caret)
library(ggplot2)


#####script to test differing # of svs in test set 


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

################
#test how different # of svs in test affect train 


missclass_by_sample<-function(id_name, test_df){
  samp<-test_df[test_df$patient_id==id_name]
  class_correct<-length(which(samp$missclassified==0))
  return(data.table(prop_class=(class_correct/nrow(samp)),
                    tumor_type=substr(samp$project_code[1], 1, nchar(samp$project_code[1])-3)))
}


train_tumor_out<-function(tt){
  print(tt)
 
  n_train<-10000
  df_train_minustt<-df_train_all[!(df_train_all$project_code==tt)]
  train <- df_train_minustt[sample(seq_len(nrow(df_train_minustt)), size = n_train), ]
  test_tt<-df_test_all[df_test_all$project_code==tt]

  #scale features in train and test set 
  train_scaled<-(train[, lapply(.SD, scale), .SDcols = features_toscale])
  train_scaled<-cbind(train_scaled, sv_class=train$sv_class)
  
  test_scaled<-(test_tt[, lapply(.SD, scale), .SDcols = features_toscale])
  test_scaled<-cbind(test_scaled, sv_class=test_tt$sv_class)
  
  
  classifier = svm(formula = sv_class ~ .,
                   data = train_scaled, type = 'C-classification',
                   kernel = 'linear', cachesize=400, 
                   probability=T)
  
  
  #predict test set with classifier 
  y_pred_all <- predict(classifier, newdata = test_scaled, decision.values = T, probability = T)
  probabilities_all<-data.table(attr(y_pred_all, 'probabilities'))
  setcolorder(probabilities_all, c('0', '1'))
  pr_all<-prediction(as.numeric(probabilities_all$`1`), as.numeric(test_scaled$sv_class))
  auc_all <- performance(pr_all, measure = "auc")
  
  
  test_tt_res<-cbind(test_tt, y_pred_all)
  split_name<-strsplit(test_tt_res$sample, '[.]')
  short_name<-lapply(1:length(split_name), function(i){
    return(split_name[[i]][1])
  })
  test_tt_res$patient_id<-short_name
  test_tt_res[,missclassified:=ifelse(sv_class==y_pred_all, 0, 1)]
  
  missclass_tt<-rbindlist(lapply(unique(test_tt_res$patient_id), missclass_by_sample, test_df=test_tt_res))
  

  return(cbind(auc_all=auc_all@y.values[[1]], missclass_tt))
  
}

project_codes<-unique(df_test_all$project_code)
tt_df<-rbindlist(lapply(project_codes[!(project_codes=='READ-US')], train_tumor_out))
#tt_df<-rbindlist(lapply(project_codes[3:5], train_tumor_out))

#write.csv(tt_df, '/Users/shu/germline_svm/data/20220615_train_minusonetumor_10ktrain.csv')
  
  
tt_df<-fread('/Users/shu/germline_svm/data/20220615_train_minusonetumor_10ktrain.csv')


tt_med<-aggregate(prop_class ~ tumor_type, tt_df, 'median')
tt_mean<-aggregate(prop_class ~ tumor_type, tt_df, 'mean')
tt_sd<-aggregate(prop_class ~ tumor_type, tt_df, 'sd')

tt_summary<-cbind(tt_med, tt_mean[,2], tt_sd[,2])

######################

#ggplot scale
scale<-1.5
theme_ss <- theme_bw(base_size=16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.text.x = element_text(angle = 90, vjust = 0.5, size=8*scale, family="mono"),
        axis.text.x = element_text(angle=45, vjust = 0.5, size=8*scale),
        #axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y = element_text(hjust = 0.5,size=10*scale),
        #axis.text = element_text(size = 10*scale, family = "mono"))
        plot.title=element_text(size=14*scale, hjust=0.5, color="black"),
        plot.subtitle=element_text(size=8*scale, hjust=0.5, face="italic", color="#252D2D")
  )



#make boxplot to look at #SVs per tumor type, and soomatic/germline proportion 

p<-ggplot(tt_df, aes(x=tumor_type, y=prop_class, fill=tumor_type)) +
  #geom_boxplot() + geom_point(aes(), size = 1, shape = 21) + ylim(c(0,1)) + 
  #scale_y_continuous(trans='log10') + 
  geom_violin() + geom_point(aes(), size = 1, shape = 21) + ylim(c(0,1)) + 
  theme_ss
p<-p + labs(title='Training on All Tumor Types Except Test Tumor Type', x='Tumor Type', y='Proportion of Correctly Classified SVs')  +  
  guides(fill=guide_legend(title="Tumor Type"))
pdf('/Users/shu/germline_svm/figs/20220615_train_minusonetumor_10ktrain_violin.pdf', width=10, height=5)
p
dev.off()







