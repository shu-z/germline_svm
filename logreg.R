library(data.table)
library(e1071)
library(caTools)
library(ROCR)
library(stringr)
library(ggfortify)
library(corrplot)



#use after train/test sets made from svm_main 

#look at correlations of features 
continuous_vars<-c('homlen', 'insertion_len', 'SPAN', 'gnomad_dist', 'hom_gc', 
                   'insertion_gc', 'line_dist', 'sine_dist', 'num_sv_sample')

somatic_features<-somatic_subset[, ..continuous_vars]
cor_somatic <- cor(somatic_features, method='pearson')
germline_features<-germline_subset[, ..continuous_vars]
cor_germline <- cor(germline_features, method='pearson')


#this to get the pvalues
library(Hmisc)
cor_somatic <-rcorr(as.matrix(somatic_features), type='spearman')
View(cor_somatic$r)
View(cor_somatic$P)
cor_germline <- rcorr(as.matrix(germline_features), type='spearman')
V#iew(cor_somatic$r)
View(cor_germline$P)
cor_all <- rcorr(as.matrix(df[, ..features_toscale]), type='spearman')
View(cor_all$P)


pdf('/Users/shu/germline_svm/figs/20220426_newannot_corrplot_pearson_continuousfeatures.pdf', width=8, height=8)
corrplot(cor_somatic, type='upper', addCoef.col = 'grey50', number.cex = 0.75, mar=c(0,0,2,0)) 
mtext("somatic, pearson, continuous features", line=2, cex=1)
corrplot(cor_germline, type='upper', addCoef.col = 'grey50', number.cex = 0.75, mar=c(0,0,2,0))
mtext("germline, pearson, continuous features", line=2, cex=1)
dev.off()

#look at all
#pdf('/Users/shu/germline_svm/figs/20220426_newannot_corrplotcombined_pearson_allfeatures.pdf', width=8, height=8)
corrplot(cor_all, type='upper', addCoef.col = 'grey50', number.cex = 0.75, mar=c(0,0,2,0)) 
mtext("10k germline, 10k somatic, pearson, all features", line=2, cex=1)
dev.off()

#logistic regression to all feauters
glm.fit <- glm(sv_class ~ ., data= train_scaled, family=binomial)
summary(glm.fit)

#look at predictions on test data
#glm.probs <- predict(glm.fit,type = "response")
glm.probs <- predict(glm.fit, newdata = test_scaled, type = "response")

#if we just want to classify as germline somatic
#glm.pred <- ifelse(glm.probs > 0.5, 1, 0)

test5<-cbind(pred_probs=data.table(glm.probs), test_scaled$sv_class)

#look at histogram of prediction probabilities
#pdf('/Users/shu/germline_svm/figs/20220419_logreg_20k_predprob_259samples_newannot.pdf', width=10)
hist(glm.probs, 10, main='Prediction Probabilities, Logistic Regression')
dev.off()


#roc/auc of logistic regression 
#pdf('/Users/shu/germline_svm/figs/20220419_logreg_20k_auc_259samples_newannot.pdf', width=10, height=6)
pr<-prediction(as.numeric(test5$pred_probs), as.numeric(test_scaled$sv_class))
auc_ROCR <- performance(pr, measure = "auc")
auc_val <- auc_ROCR@y.values[[1]]
pref <- performance(pr, "tpr", "fpr")
#plot(pref, main = paste0("ROC Curve", '\n', 'AUC: ', substr(auc_val, 1,6) ), colorize = F)
plot(pref, main = paste0('AUC: ', substring(auc_val,1,5), '\n', 'Logistic Regression, 16k train, 4k test'), colorize = F, 
     cex.lab=1, cex.main=1)
abline(a = 0, b = 1)
dev.off()




library(broom)
library(dplyr)


run_single_logreg<-function(i, data){
  train_df<-cbind(train_scaled[, ..i], sv_class=train_scaled[, sv_class])
  test_df<-cbind(test_scaled[, ..i], sv_class=test_scaled[, sv_class])
  
  
  glm.fit <- glm(sv_class ~ ., data= train_df, family=binomial)

  #look at predictions on test data
  glm.probs <- predict(glm.fit, newdata = test_df, type = "response")

  #roc/auc of logistic regression
  pr<-prediction(as.numeric(glm.probs), as.numeric(test_df$sv_class))
  auc <- performance(pr, measure = "auc")

  coefficient_df<-summary(glm.fit)$coefficients
  
  #plot auc all on same plot
  
  
  return(cbind(feature=colnames(train_df)[1], auc=auc@y.values[[1]], data.table(t(coefficient_df[2,]))))



}



single_logreg_df<-rbindlist(lapply(1:16, run_single_logreg, data=train_scaled))



####same as above, but for loop for sake of plotting

#pdf('/Users/shu/germline_svm/figs/singlefeature_logreg_auc.pdf', width=8, height=7)
pdf('/Users/shu/germline_svm/figs/singlefeature_logreg_predprob_hist.pdf', width=8, height=7)

colors<-cols25(16)
for (i in 1:15){
  train_df<-cbind(train_scaled[, ..i], sv_class=train_scaled[, sv_class])
  test_df<-cbind(test_scaled[, ..i], sv_class=test_scaled[, sv_class])
  
  
  glm.fit <- glm(sv_class ~ ., data= train_df, family=binomial)
  
  #look at predictions on test data
  glm.probs <- predict(glm.fit, newdata = test_df, type = "response")
  coefficient_df<-summary(glm.fit)$coefficients
  hist(glm.probs, seq(0,1,length.out=100), main=paste0(features_toscale[i]), xlim=c(0,1))
  
  #roc/auc of logistic regression
  pr<-prediction(as.numeric(glm.probs), as.numeric(test_df$sv_class))
  auc <- performance(pr, measure = "auc")

  print(auc@y.values[[1]])
  #plot auc all on same plot
  # pref <- performance(pr, "tpr", "fpr")
  # plot(pref, main = 'Single Feature Logistic Regression AUCs', colorize = F, 
  #      cex.lab=1, cex.main=1, add=(i!=1), lwd=2, col=colors[i])
  # abline(a = 0, b = 1)

  
  #return(cbind(feature=colnames(train_df)[1], auc=auc@y.values[[1]], data.table(t(coefficient_df[2,]))))
}

#legend(0.79, 0.52, legend= features_toscale, fill=colors,  cex=0.8)
dev.off()


#########

####try single svm see if anything different 
run_single_svm<-function(i, data){
  train_df<-cbind(train_scaled[, ..i], sv_class=train_scaled[, sv_class])
  test_df<-cbind(test_scaled[, ..i], sv_class=test_scaled[, sv_class])
  
  classifier = svm(formula = sv_class ~ .,
                   data = train_df, type = 'C-classification',
                   kernel = 'radial', cachesize=400, 
                   probability=T)
  
  #predict test set with classifier 
  y_pred <- predict(classifier, newdata = test_df, decision.values = T, probability = T)
  probabilities<-attr(y_pred, 'probabilities')
  
  max_prob<-apply(probabilities,1, function(x){
    val<-as.vector(x)
    return(val[2])
  })
  
  pr<-prediction(as.numeric(max_prob), as.numeric(test_df$sv_class))
  auc <- performance(pr, measure = "auc")

  return(data.table(feature=colnames(train_df)[1], auc=auc@y.values[[1]]))
  
  
  
}

single_svm_df<-rbindlist(lapply(1:16, run_single_svm, data=train_scaled))


write.csv(single_svm_df, '/Users/shu/germline_svm/data/20220504_single_feature_svm_16ktrain.csv')
write.csv(single_logreg_df, '/Users/shu/germline_svm/data/20220504_single_feature_logreg_16ktrain.csv')



