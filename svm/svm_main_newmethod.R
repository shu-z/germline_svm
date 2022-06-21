library(data.table)
library(e1071)
library(caTools)
library(ROCR)
library(stringr)
library(ggfortify)
library(caret)
library(ggplot2)


#use after add_features, once somatic and germline subsets are made

##########
#make train and test sets

#df_train_all is all somatic SVs from train set, with equal number of germline SVs
#with 20220601 sample, 62,000 SVs total 
somatic_train<-somatic_subset[somatic_subset$train_test=='train']
germline_train<-germline_subset[germline_subset$train_test=='train']
df_train_all<-rbind(somatic_train, germline_train[sample(seq_len(nrow(germline_train)), size = nrow(somatic_train)),])


features_toscale<-c('homlen', 'insertion_len', 'SPAN', 'gnomad_dist', 'del', 'dup', 'inv', 'inter',
                    #'hom_gc', 'insertion_gc', 'line_dist', 'sine_dist', 'num_sv_sample', 'CN_annot', 'exon_annot')
                     'line_dist', 'sine_dist', 'num_sv_sample', 'CN_annot', 'exon_annot')


#df_test_all is all long SVs from test tumors (10 per tumor type)
somatic_test<-somatic_subset[somatic_subset$train_test=='test']
germline_test<-germline_subset[germline_subset$train_test=='test']
df_test_all<-rbind(somatic_test, germline_test)


#if we want our test set to have ~equal SVs per tumor 
even_sv_sampling<-function(sample){
  long_df<-df_test_all[df_test_all$sample_name==sample_name]
  n_len<-min(nrow(long_df), 600)
  return((long_df[sample(seq_len(nrow(long_df)), size = n_len),]))
}
df_test_even<-rbindlist(lapply(unique(df_test_all$sample_name), even_sv_sampling))



################
n_train<-30000

#define train and test sets
train <- df_train_all[sample(seq_len(nrow(df_train_all)), size = n_train), ]
#test <- df_test_all[sample(seq_len(nrow(df_test_all)), size = n_test), ]
test<-df_test_all

#scale features in train and test set 
train_scaled<-(train[, lapply(.SD, scale), .SDcols = features_toscale])
train_scaled<-cbind(train_scaled, sv_class=train$sv_class)
test_scaled<-(test[, lapply(.SD, scale), .SDcols = features_toscale])
test_scaled<-cbind(test_scaled, sv_class=test$sv_class)

#run classifier (linear kernel)
classifier = svm(formula = sv_class ~ .,
                 data = train_scaled, type = 'C-classification',
                 kernel = 'linear', cachesize=400, 
                 probability=T)


# classifier_radial = svm(formula = sv_class ~ .,
#                  data = train_scaled, type = 'C-classification',
#                  kernel = 'radial', cachesize=400, 
#                  probability=T)


#tune classifier (select gamma and C ranges)
tune_obj <- tune(svm, sv_class~., data = train_scaled,
            #ranges = list(gamma = 2^(-6:2), cost = 2^(0:4)),
            ranges = list(gamma = 2^(-1:2), cost = 10^(1:5)),
            tunecontrol = tune.control(sampling = "fix"))
summary(tune_obj)

#save tuning output
#need to rerun classifier with best parameters
#saveRDS(tune_obj, '/Volumes/xchip_beroukhimlab/Shu/ccle/poster_figs/svm_tuning_obj.rds')
#pdf('/Volumes/xchip_beroukhimlab/Shu/ccle/poster_figs/svm_performance.pdf')
plot(tune_obj)
#dev.off()


#predict test set SVs with classifier
y_pred <- predict(classifier, newdata = test_scaled, decision.values = T, probability = T)
probabilities<-data.table(attr(y_pred, 'probabilities'))
setcolorder(probabilities, c('0', '1'))

#get performance of classifier
pr<-prediction(as.numeric(probabilities$`1`), as.numeric(test_scaled$sv_class))
#pr<-prediction(as.numeric(y_pred), as.numeric(test_scaled$sv_class))
auc <- performance(pr, measure = "auc")


#graph roc/auc
#pdf('/Volumes/xchip_beroukhimlab/Shu/ccle/poster_figs/20211214_roc_50k_10features_wide.pdf', width=10, height=6)
auc_val <- auc@y.values[[1]]
pref <- performance(pr, "tpr", "fpr")
#plot(pref, main = paste0("ROC curve", '\n', 'AUC: ', auc_val), colorize = F)
plot(pref, main = paste0('AUC: ', substring(auc_val,1,5) ), colorize = F, 
     cex.lab=1.5, cex.main=1.5)
abline(a = 0, b = 1)
#dev.off()
  
  
  
#platt's calibration for prediction probabilities
source('/Users/shu/GermlineSVAnnotator/svm/prcalibrate_eric.R')

#pdf('/Users/shu/germline_svm/figs/202204019_svm_20k_259samples_newannot.pdf', width=5)
test_calibrate<-prCalibrate(as.numeric(levels(test_scaled$sv_class))[test_scaled$sv_class], as.numeric(max_prob))
pr_cal<-prediction((test_calibrate$cal.probs), as.numeric(test_scaled$sv_class))

#plot new auc curve
auc_ROCR_cal <- performance(pr_cal, measure = "auc")
auc_val_cal <- auc_ROCR_cal@y.values[[1]]
pref_cal <- performance(pr_cal, "tpr", "fpr")
#plot(pref, main = paste0("ROC curve", '\n', 'AUC: ', auc_val), colorize = F)
plot(pref_cal, main = paste0('AUC: ', substring(auc_val_cal,1,5) ), colorize = F, 
     cex.lab=1.5, cex.main=1.5)
abline(a = 0, b = 1)
dev.off()
  
  
  
######logistic regression stuff to compare to svm

#run logistic regression, using scaled train set
glm.fit <- glm(sv_class ~ ., data= train_scaled, family=binomial)
summary(glm.fit)

#look at predictions on test data
glm.probs <- predict(glm.fit, newdata = test_scaled, type = "response")

#if we just want to classify as germline somatic
glm.pred <- ifelse(glm.probs > 0.5, 1, 0)

#bind prediction probabilities to test dt
logreg_test<-cbind(pred_probs=data.table(glm.probs), test_scaled$sv_class)

#look at histogram of prediction probabilities
#pdf('/Users/shu/germline_svm/figs/20220419_logreg_20k_predprob_259samples_newannot.pdf', width=10)
hist(glm.probs, 30, main='Prediction Probabilities, Logistic Regression')
#dev.off()


#roc/auc of logistic regression 
#pdf('/Users/shu/germline_svm/figs/20220419_logreg_20k_auc_259samples_newannot.pdf', width=10, height=6)
pr_logreg<-prediction(as.numeric(logreg_test$pred_probs), as.numeric(test_scaled$sv_class))
auc_logreg <- performance(pr_logreg, measure = "auc")
auc_val_logreg <- auc_logreg@y.values[[1]]
print(auc_val_logreg)
pref_logreg <- performance(pr_logreg, "tpr", "fpr")
#plot(pref, main = paste0("ROC Curve", '\n', 'AUC: ', substr(auc_val, 1,6) ), colorize = F)
plot(pref_logreg, main = paste0('AUC: ', substring(auc_val_logreg,1,5), '\n', 'Logistic Regression'), colorize = F, 
     cex.lab=1, cex.main=1)
abline(a = 0, b = 1)
dev.off()
  
  





