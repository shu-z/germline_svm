library(data.table)
library(e1071)
library(caTools)
library(ROCR)
library(stringr)
library(ggfortify)
library(caret)
library(ggplot2)
library(scales)




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
#get probabilities 
y_pred <- predict(classifier, newdata = test_scaled, decision.values = T, probability = T)
probabilities<-data.table(attr(y_pred, 'probabilities'))
setcolorder(probabilities, c('0', '1'))
test_withprob<-cbind(test, prob=probabilities$`1`)



test_pred_cutoffs<-function(cutoff){
  
  test_predprob<-test_withprob[,c('sv_class', 'prob')]
  test_predprob[, pred:=ifelse(prob<cutoff, 0, 1)]
  test_predprob[, misclass:=ifelse(sv_class==pred, 0,1)]
  
  #split test set by germline somatic 
  test_germline<-test_predprob[test_predprob$sv_class==0]
  test_somatic<-test_predprob[test_predprob$sv_class==1]
  
  #calculate values 
  prop_class<-sum(test_predprob$misclass==0)/nrow(test_predprob)
  
  #calculate proportion of germline/somatic called correctly  
  prop_germline_class<-1-(sum(test_germline$misclass==1)/nrow(test_germline))
  prop_somatic_class<-1-(sum(test_somatic$misclass==1)/nrow(test_somatic))
    
  #calculate true proportion of germline/somatic SVs in called germline/somatic SVs
  prop_true_germline<-(sum(test_germline$misclass==0))/sum(test_predprob$pred==0)
  prop_true_somatic<-(sum(test_somatic$misclass==0))/sum(test_predprob$pred==1)
  

  return(data.table(cutoff, prop_class, prop_germline_class, prop_somatic_class, 
                    prop_true_germline, prop_true_somatic))
  
  
}




cutoff_df<-rbindlist(lapply(seq(0,1, by=0.1), test_pred_cutoffs))

###################plots 

#plot scale
scale<-1.5
theme_ss <- theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(vjust = 0.5, size=8*scale),
        #axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y = element_text(hjust = 0.5,size=10*scale)
        #axis.text = element_text(size = 10*scale, family = "mono"))
  )




#plot proportion of SVs classified correctly 
cutoff_df_plot<-melt(cutoff_df[,c(1:4)], id='cutoff', variable.name = 'SV_Class')

p<-ggplot(cutoff_df_plot, aes(x=cutoff, y=value, group=SV_Class)) + 
  geom_line(aes(color=SV_Class)) + geom_point(aes(color=SV_Class), size=0.5) 
p<- p + xlab('Prediction Probability Cutoff') + ylab ('Proportion of Misclassified SVs') + 
  scale_colour_manual(labels=c('All SVs', 'Germline SVs', 'Somatic SVs'), 
                      values=c( '#6D72C3', '#F8766D', '#00BFC4'))  + 
 scale_x_continuous(breaks = seq(0,1, by=0.2)) + theme_ss  

pdf('/Users/shu/germline_svm/figs/20220622_svmlinear_30k_predprobcutoff.pdf', width=8, height=5)
p
dev.off()





#plot true positives 
truepositives_df<-melt(cutoff_df[,c(1,5,6)], id='cutoff', variable.name = 'SV_Class')
p<-ggplot(truepositives_df, aes(x=cutoff, y=value, group=SV_Class)) + 
  geom_line(aes(color=SV_Class)) + geom_point(aes(color=SV_Class), size=0.5) 
p<- p + xlab('Prediction Probability Cutoff') + ylab ('True SVs/Classified SVs') + 
  scale_colour_manual(labels=c('Germline SVs', 'Somatic SVs'),
                     values=c( '#F8766D', '#00BFC4'))  +
  scale_x_continuous(breaks = seq(0,1, by=0.2)) + theme_ss  

pdf('/Users/shu/germline_svm/figs/20220622_svmlinear_30k_predprobcutoff_truepositives.pdf', width=8, height=5)
p
dev.off()



 








