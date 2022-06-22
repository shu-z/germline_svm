library(data.table)
library(e1071)
library(caTools)
library(ROCR)
library(stringr)
library(ggfortify)
library(caret)
library(ggplot2)


#####script to test differing # 


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


####################
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


################
#test how changing prediction cutoff affects number of germline and somatic SVs identified


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




