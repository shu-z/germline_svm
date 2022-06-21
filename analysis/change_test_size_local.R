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


# features_toscale<-c('homlen', 'insertion_len', 'SPAN', 'gnomad_dist', 'del', 'dup', 'inv', 'inter',
#                     'hom_gc', 'insertion_gc', 'line_dist', 'sine_dist', 'num_sv_sample', 'CN_annot', 'exon_annot')

features_toscale<-c('homlen', 'insertion_len', 'SPAN', 'gnomad_dist', 'del', 'dup', 'inv', 'inter',
                     'line_dist', 'sine_dist', 'num_sv_sample', 'CN_annot', 'exon_annot')

somatic_test<-somatic_subset[somatic_subset$train_test=='test']
germline_test<-germline_subset[germline_subset$train_test=='test']
df_test_all<-rbind(somatic_test, germline_test)


#if we want our test set to have ~equal SVs per tumor 
even_sv_sampling<-function(samp_name){
  long_df<-df_test_all[df_test_all$sample_name==samp_name]
  n_len<-min(nrow(long_df), 600)
  return((long_df[sample(seq_len(nrow(long_df)), size = n_len),]))
}


df_test_evensamp<-rbindlist(lapply(unique(df_test_all$sample_name), even_sv_sampling))
df_test_somatic<-df_test_evensamp[df_test_evensamp$CLASS=='SOMATIC']
df_test_germline<-df_test_evensamp[df_test_evensamp$CLASS=='GERMLINE']

df_test_samp1<-rbind(df_test_somatic[sample(seq_len(nrow(df_test_somatic)), size = 10000),],
                     df_test_germline[sample(seq_len(nrow(df_test_germline)), size = 10000),])

################
#test how different # of svs in test affect train 

change_n_train_svm<-function(n_train){
  print(n_train)
  #test on 1000 SVs for now
  #n_test<-10000
  #next sample for testing, different SV proportions 
  train <- df_train_all[sample(seq_len(nrow(df_train_all)), size = n_train), ]
  test_all<-df_test_all
  test_evensamp<-df_test_samp1
  
  #scale features in train and test set 
  train_scaled<-(train[, lapply(.SD, scale), .SDcols = features_toscale])
  train_scaled<-cbind(train_scaled, sv_class=train$sv_class)

  
  test_scaled<-(test_all[, lapply(.SD, scale), .SDcols = features_toscale])
  test_scaled<-cbind(test_scaled, sv_class=test_all$sv_class)
  
  test_even_scaled<-(test_evensamp[, lapply(.SD, scale), .SDcols = features_toscale])
  test_even_scaled<-cbind(test_even_scaled, sv_class=test_evensamp$sv_class)
  
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
  
  
  #repeat for evenly sampled svs
  y_pred_even <- predict(classifier, newdata = test_even_scaled, decision.values = T, probability = T)
  probabilities_even<-data.table(attr(y_pred_even, 'probabilities'))
  setcolorder(probabilities_even, c('0', '1'))
  
  pr_even<-prediction(as.numeric(probabilities_even$`1`), as.numeric(test_even_scaled$sv_class))
  auc_even <- performance(pr_even, measure = "auc")
  
  
  
  return(data.table(n_train, auc_all=auc_all@y.values[[1]], auc_even_sv=auc_even@y.values[[1]]))
  
}


for (i in 1:10){
  print(paste0('i: ', i))

  n_train_df<-rbindlist(lapply(c(100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000, 20000, 30000, 40000, 50000, 60000), change_n_train_svm))
  n_train_df<-cbind(n_train_df, i=i)
  
  #write.csv(n_train_df, paste0('/Users/shu/germline_svm/data/change_train_size/13_feat/20220609_svmlinear_change_trainsize_i', i, '.csv'))
  
  
}



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



#now, read in all dfs and combine

all_training_list<-list()
for (i in 1:10){
  #one_iter<-fread(paste0('/Users/shu/germline_svm/data/change_train_size/13_feat/20220601_change_trainsize_i', i, '.csv'), drop=1)
  #one_iter<-fread(paste0('/Users/shu/germline_svm/data/change_train_size/all_feat/20220606_change_trainsize_i', i, '.csv'), drop=1)
  
  #this one for log reg 
  #one_iter<-fread(paste0('/Users/shu/germline_svm/data/change_train_size/13_feat/20220609_logreg_change_trainsize_i', i, '.csv'), drop=1)
  #one_iter<-fread(paste0('/Users/shu/germline_svm/data/change_train_size/all_feat/20220606_logreg_change_trainsize_i', i, '.csv'), drop=1)
  
  #svm_linear
  one_iter<-fread(paste0('/Users/shu/germline_svm/data/change_train_size/13_feat/20220609_svmlinear_change_trainsize_i', i, '.csv'), drop=1)
  #one_iter<-fread(paste0('/Users/shu/germline_svm/data/change_train_size/all_feat/20220609_svmlinear_change_trainsize_i', i, '.csv'), drop=1)
  
  
  
  
  all_training_list[[i]]<-one_iter
}

all_training_df<-rbindlist(all_training_list)
all_training_df[,n_train:=as.factor(n_train)]
all_training_df[,i:=as.factor(i)]

#get mean values across training set size
med_col<-aggregate(auc_all ~ n_train, all_training_df, 'median')
mean_col<-aggregate(auc_all ~ n_train, all_training_df, 'mean')
sd_col<-aggregate(auc_all ~ n_train, all_training_df, 'sd')
all_df_stats<-data.table(cbind(med_col, mean_col[,2], sd_col[,2]))
colnames(all_df_stats)<-c('n', 'median', 'mean', 'sd')
all_df_stats[,coef_var:= (sd/mean)]


#make boxplot to look at #SVs per tumor type, and soomatic/germline proportion 

p<-ggplot(all_training_df, aes(x=n_train, y=auc_all, fill=i)) +
  geom_point(aes(), size = 1, shape = 21) + 
  stat_summary(aes(y = auc_all,group=1), fun=base::mean, colour="grey", geom="line") + 
  stat_summary(aes(y = auc_all, group=1, width=0.3),fun.data= mean_se, colour='grey', geom = "errorbar") + 
  #scale_y_continuous(trans='log10') + 
  #ylim(c(0.83,0.95)) + theme_ss
  ylim(c(0.8,1)) + theme_ss

  #ylim(c(0.4,1)) + theme_ss
  
p<-p + labs(title='Training Set Size vs AUC (Radial SVM)', x='Training Set Size', y='AUC',
            subtitle=paste0('330 samples in training set, 200 samples in test set', '\n', 
            #'Test on 10k germline, 10k somatic non-short SVs', '\n', 
            'Test on all non-short SVs', '\n',
            #'Homology GC and Insertion GC removed'))  +  
            'All Features'))  +  
            
  guides(fill=guide_legend(title="Iteration"))
#pdf('/Users/shu/germline_svm/figs/change_train_size/20220610_svmradial_changetrain_allfeat_zoom.pdf', width=10, height=7)
p
dev.off()







####


##aggregate data 

dir_path<-'/Users/shu/germline_svm/data/change_train_size/'
paths<-c('13_feat/20220601_change_trainsize_i', 'all_feat/20220606_change_trainsize_i',
         '13_feat/20220609_svmlinear_change_trainsize_i', 'all_feat/20220609_svmlinear_change_trainsize_i',
         '13_feat/20220609_logreg_change_trainsize_i','all_feat/20220606_logreg_change_trainsize_i')

aggregate_data<-function(path){
  all_training_list<-list()
  for (i in 1:10){
    one_iter<-fread(paste0(dir_path, paths[path], i, '.csv'), drop=1)
    all_training_list[[i]]<-one_iter
  }
    
    all_training_df<-rbindlist(all_training_list)
    all_training_df[,n_train:=as.factor(n_train)]
    all_training_df[,i:=as.factor(i)]
    
    #get mean values across training set size
    mean<-aggregate(auc_all ~ n_train, all_training_df, 'mean')
    median<-aggregate(auc_all ~ n_train, all_training_df, 'median')
    sd<-aggregate(auc_all ~ n_train, all_training_df, 'sd')
  
    df<-cbind(mean, median[,2], sd[,2])
    colnames(df)<-c('i', 'auc_mean', 'auc_med', 'auc_sd')
    return(df)


}


test<-(lapply(1:6, aggregate_data))

df_agg<-(test[[1]])
for (i in c(3,5)){
  df_i<-test[[i]]
  df_agg<-cbind(df_agg, df_i[,2:4])
}

colnames(df_agg)<-c('i', 'svm_radial_mean', 'svm_radial_median', 'svm_radial_sd',
                 'svm_linear_mean', 'svm_linear_median', 'svm_linear_sd', 
                 'logreg_mean', 'logreg_median', 'logreg_sd')

df_agg<-data.table(df_agg)
df_agg[, linear_coefvar:=svm_linear_sd/svm_linear_mean]

#write.csv(df_agg, '/Users/shu/germline_svm/data/change_train_size/allfeat_summary_stats.csv')

###################### do for logreg



change_n_train_logreg<-function(n_train){
  print(n_train)
  #test on 1000 SVs for now
  #n_test<-10000
  #next sample for testing, different SV proportions 
  train <- df_train_all[sample(seq_len(nrow(df_train_all)), size = n_train), ]
  test_all<-df_test_all
  test_evensamp<-df_test_samp1
  
  #scale features in train and test set 
  train_scaled<-(train[, lapply(.SD, scale), .SDcols = features_toscale])
  train_scaled<-cbind(train_scaled, sv_class=train$sv_class)
  
  
  test_scaled<-(test_all[, lapply(.SD, scale), .SDcols = features_toscale])
  test_scaled<-cbind(test_scaled, sv_class=test_all$sv_class)
  
  test_even_scaled<-(test_evensamp[, lapply(.SD, scale), .SDcols = features_toscale])
  test_even_scaled<-cbind(test_even_scaled, sv_class=test_evensamp$sv_class)
  
  
  
  glm.fit <- glm(sv_class ~ ., data= train_scaled, family=binomial)
  summary(glm.fit)
  
  #if we just want to classify as germline somatic
  #glm.pred <- ifelse(glm.probs > 0.5, 1, 0)
  
  
  #predict test set with classifier 
  y_pred_all <- predict(glm.fit, newdata = test_scaled, type = "response")
  pr_all<-prediction(as.numeric(y_pred_all), as.numeric(test_scaled$sv_class))
  auc_all <- performance(pr_all, measure = "auc")
  
  
  #repeat for evenly sampled svs
  y_pred_even <- predict(glm.fit, newdata = test_even_scaled, type = "response")
  pr_even<-prediction(as.numeric(y_pred_even), as.numeric(test_even_scaled$sv_class))
  auc_even <- performance(pr_even, measure = "auc")
  
  
  
  return(data.table(n_train, auc_all=auc_all@y.values[[1]], auc_even_sv=auc_even@y.values[[1]]))
  
}


for (i in 1:10){
  print(paste0('i: ', i))
  n_train_df<-rbindlist(lapply(c(100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000, 20000, 30000, 40000, 50000, 60000), change_n_train_logreg))
  #n_train_df<-rbindlist(lapply(c(100, 250, 500, 750), change_n_train_svm))
  
  n_train_df<-cbind(n_train_df, i=i)
  write.csv(n_train_df, paste0('/Users/shu/germline_svm/data/change_train_size/13_feat/20220609_logreg_change_trainsize_i', i, '.csv'))
  
  
}

