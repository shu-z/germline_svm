
library(ggplot2)


#TEST DIFFERENT PROPORTIONS OF GERMLINE/SOMATIC IN TEST SET YIPPEPEPEPEPEP
#for use after features are added 



#to test different proportions of germline/somatic in train and test sets 
somatic_subset<-df[df$CLASS=='SOMATIC']
germline_subset<-df[df$CLASS=='GERMLINE']

#first sample for training 
n_train<-10000
train_somatic_ind <- sample(seq_len(nrow(somatic_subset)), size = n_train)
train_germline_ind <- sample(seq_len(nrow(germline_subset)), size = n_train)
train<-rbind(somatic_subset[train_somatic_ind,], germline_subset[train_germline_ind,])

# feature scaling
features_toscale<-c('homlen', 'insertion_len', 'SPAN', 'gnomad_dist', 'del', 'dup', 'h2hINV', 't2tINV', 'inter',
                    'hom_gc', 'insertion_gc', 'line_dist', 'sine_dist', 'num_sv_sample', 'CN_annot', 'exon_annot')
train_scaled<-(train[, lapply(.SD, scale), .SDcols = features_toscale])
train_scaled<-cbind(train_scaled, sv_class=train$sv_class)

#run classifier
classifier = svm(formula = sv_class ~ .,
                 data = train_scaled, type = 'C-classification',
                 kernel = 'radial', cachesize=400, 
                 probability=T)


#remove SVs that became part of train set 
somatic_subset <- somatic_subset[-train_somatic_ind, ]
germline_subset <- germline_subset[-train_germline_ind, ]



change_test_prop<-function(prop_somatic, data){
  
  #next sample for testing, different SV proportions 
  #use 10k for total test set for now 
  n_test<-10000
  test_somatic_ind <- sample(seq_len(nrow(somatic_subset)), size = floor(prop_somatic * n_test))
  test_germline_ind <- sample(seq_len(nrow(germline_subset)), size = floor((1-prop_somatic) * n_test))
  test<-rbind(somatic_subset[test_somatic_ind,], germline_subset[test_germline_ind,])
  test_scaled<-(test[, lapply(.SD, scale), .SDcols = features_toscale])
  test_scaled<-cbind(test_scaled, sv_class=test$sv_class)
  
  
  #predict test set with classifier 
  y_pred <- predict(classifier, newdata = test_scaled, decision.values = T, probability = T)
  probabilities<-attr(y_pred, 'probabilities')
  
  max_prob<-apply(probabilities,1, function(x){
    val<-as.vector(x)
    return(val[1])
  })
  
  pr<-prediction(as.numeric(max_prob), as.numeric(test_scaled$sv_class))
  auc <- performance(pr, measure = "auc")

  return(data.table(prop_somatic, auc=auc@y.values[[1]]))
  
}

prop_class_svm_df<-rbindlist(lapply(seq(0.1, 0.9, by=0.1), change_test_prop))




#TO TEST BOTH DIFFERENCES IN TRAIN AND TEST

features_toscale<-c('homlen', 'insertion_len', 'SPAN', 'gnomad_dist', 'del', 'dup', 'h2hINV', 't2tINV', 'inter',
                    'hom_gc', 'insertion_gc', 'line_dist', 'sine_dist', 'num_sv_sample', 'CN_annot', 'exon_annot')

#test function, to be called by each train function 
change_test_prop<-function(prop_somatic_test, somatic_subset_test, germline_subset_test, classifier){
  
  #next sample for testing, different SV proportions 
  #use 10k for total test set for now 
  n_test<-10000
  test_somatic_ind <- sample(seq_len(nrow(somatic_subset_test)), size = floor(prop_somatic_test * n_test))
  test_germline_ind <- sample(seq_len(nrow(germline_subset_test)), size = floor((1-prop_somatic_test) * n_test))
  test<-rbind(somatic_subset_test[test_somatic_ind,], germline_subset_test[test_germline_ind,])
  test_scaled<-(test[, lapply(.SD, scale), .SDcols = features_toscale])
  test_scaled<-cbind(test_scaled, sv_class=test$sv_class)
  
  
  #predict test set with classifier 
  print(dim(test_scaled))
  y_pred <- predict(classifier, newdata = test_scaled, decision.values = T, probability = T)
  probabilities<-attr(y_pred, 'probabilities')
  
  max_prob<-apply(probabilities,1, function(x){
    val<-as.vector(x)
    return(val[1])
  })
  
  pr<-prediction(as.numeric(max_prob), as.numeric(test_scaled$sv_class))
  auc <- performance(pr, measure = "auc")
  
  return(data.table(prop_somatic_test, auc=auc@y.values[[1]]))
  
}


change_train_prop<-function(prop_somatic_train, somatic_subset, germline_subset){
  print(paste0('train: ', prop_somatic_train))
  #first sample for training 
  n_train<-10000
  train_somatic_ind <- sample(seq_len(nrow(somatic_subset)), size = floor(n_train*prop_somatic_train))
  train_germline_ind <- sample(seq_len(nrow(germline_subset)), size = floor(n_train*(1-prop_somatic_train)))
  train<-rbind(somatic_subset[train_somatic_ind,], germline_subset[train_germline_ind,])
  
  # feature scaling
  train_scaled<-(train[, lapply(.SD, scale), .SDcols = features_toscale])
  train_scaled<-cbind(train_scaled, sv_class=train$sv_class)
  #run classifier
  classifier = svm(formula = sv_class ~ .,
                   data = train_scaled, type = 'C-classification',
                   kernel = 'radial', cachesize=400, 
                   probability=T)
  
  
  #remove SVs that became part of train set 
  somatic_subset_test <- somatic_subset[-train_somatic_ind, ]
  germline_subset_test <- germline_subset[-train_germline_ind, ]
  
  #call test function 
  prop_class_svm_df<-rbindlist(lapply(seq(0.1, 0.9, by=0.1), change_test_prop, somatic_subset_test, germline_subset_test, classifier))
  return(cbind(prop_somatic_train, prop_class_svm_df))
}



svm_test_df<-rbindlist(lapply(seq(0.1, 0.9, by=0.1), change_train_prop, somatic_subset, germline_subset))
write.csv(svm_test_df,'/Users/shu/germline_svm/data/20220504_test_propgermline_prop_somatic_n10kboth.csv')






pdf('/Users/shu/germline_svm/figs/20220504_test_propgermline_prop_somatic_n10kboth.pdf', width=15)
ggplot(svm_test_df, aes(x=prop_somatic_test, y=auc, 
  group=prop_somatic_train, colour=factor(prop_somatic_train))) + 
  geom_line(size=0)  +  geom_smooth(se=F)+ #geom_point()
  xlab('Somatic Proportion in Test Set') +ylab("AUC") + 
  labs(color='Somatic Proportion in Train Set') + theme_classic()
  #geom_line(data = df26means, aes(group = 1), size = 1.25, color = "black")
dev.off()

library(data.table)
test<-fread('/Users/shu/germline_svm/data/20220504_test_propgermline_prop_somatic_n10kboth.csv')

aggregate(test$auc, list(test$prop_somatic_train), mean)
aggregate(test$auc, list(test$prop_somatic_train), max)





###########test distance to germline
test_set_all<-df1



all_germline<-fread('/Users/Shu/germline_svm/data/20220511_snowman_allsamples_newannot_germlineonly.csv')
all_somatic<-fread('/Users/Shu/germline_svm/data/20220511_snowman_allsamples_newannot_somaticonly.csv')

test_set_all<-rbind(all_germline, all_somatic)


test1<-test_set_all[sample(1:nrow(test_set_all), 25000),]

get_performance<-function(i, test_set){
  test_set[,sv_class:=ifelse(CLASS=='GERMLINE',0,1)]
  test_set[,sv_pred:=ifelse(gnomad_dist<i,0,1)]
  test_set<-na.omit(test_set)
  
  pr<-prediction(as.numeric(test_set$sv_pred), as.numeric(test_set$sv_class))
  tpr_fpr <- performance(pr, measure='tpr','fpr')
  tpr<-tpr_fpr@y.values[[1]][2]
  fpr<-tpr_fpr@x.values[[1]][2]
  
  
  return(data.table(tpr, fpr))
  
  
}



test1<-lapply(c(5, 10, seq(0,100, 25),seq (200,900,100), seq(1000,10000, 1000), seq(15000, 30000, 5000)), get_performance, test_set_all)
#test1<-lapply(c(seq(0,100, 25)), get_performance, test_set_all)

test2<-list(1,1)
test2<-rbind(test2, rbindlist(test1))
test2<-rbind(test2, list(0, 0))

require(pracma)
AUC = trapz(test2$fpr,test2$tpr)

pdf('/Users/shu/germline_svm/figs/20220517_distance_Roc_10k.pdf')
plot(test2$fpr, test2$tpr, main = paste0("ROC Curve from Germline Reference Distance (0-10000bp)", '\n', 'AUC: ', substr(-AUC, 1,5) ), 
     xlim=c(0,1), ylim=c(0,1), xlab='False Positive Rate', ylab='True Positive Rate', colorize = F)
lines(test2$fpr, test2$tpr)
abline(a = 0, b = 1)

dev.off()

