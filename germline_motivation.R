library(data.table)
library(e1071)
library(caTools)
library(ROCR)
library(stringr)
library(ggfortify)
library(ggplot2)


#read in snowman germline somatic calls 
#new one has line-sine annotation
#snowman_germdist_dir<-'/Volumes/xchip_beroukhimlab/Alex/CCLE_PCAWG/Germline_Correctiion/Test_Set/Germ_dist/'
snowman_germdist_dir<-'/Volumes/xchip_beroukhimlab/Alex/CCLE_PCAWG/Germline_Correctiion/Test_Set/Germ_dist/FUNCT_ANNOT/Line_Sine_Annot/'
all_paths<-paste0(snowman_germdist_dir, list.files(snowman_germdist_dir))
subset_snowman_germdist<-rbindlist(lapply(all_paths[1:100], function(file_name){
  return(cbind(file_name, fread(file_name)))}))



#source()

#df<-subset_snowman_germdist
df1<-data.table(subset_snowman_germdist)




remove_filter<-function(cutoff,df_i){
  df_cutoff<-df_i[df_i$germline_dist<cutoff]
  prop_germline<-sum(df_cutoff$CLASS=='GERMLINE')/sum(df_i$CLASS=='GERMLINE')
  prop_somatic<-sum(df_cutoff$CLASS=='SOMATIC')/sum(df_i$CLASS=='SOMATIC')
  #return(data.table(prop=c(prop_germline, prop_somatic), sv_class=c('germline', 'somatic'), i=c(i,i), sample=c(df$sample, df$sample)))
  return(data.table(prop_germline, prop_somatic,cutoff))
  
}

filter_sample<-function(x, all_df){
  df2<-all_df[all_df$file_name==x]
  results<-rbindlist(lapply(c(50,250,1000,5000), remove_filter, df_i=df2))
  return(results)
}

test<-rbindlist(lapply(unique(df1$file_name), filter_sample, all_df=df1))


calc_means<-function(x,df){
  df_i<-df[df$cutoff==x]
  print(nrow(df_i))
  return(data.table(i=c(x,x), prop_mean=c(mean(df_i$prop_germline), mean(df_i$prop_somatic)), 
                    class=c('Germline', 'Somatic'), sd=c(sd(df_i$prop_germline),sd(df_i$prop_somatic))))
  
}

test1<-rbindlist(lapply(c(50,250,1000,5000), calc_means, test))

scale<-1
theme_ss <- theme_classic() +
  #theme(axis.text.x = element_text(angle = 45, vjust = -0.5, size=8*scale),
  theme(axis.text.x = element_text( size=20*scale),
        axis.text.y = element_text( size=20*scale),
        axis.title=element_text(size=26*scale),
        legend.title=element_text(size=30*scale),
        legend.text=element_text(size=26*scale)) 


pdf('/Volumes/xchip_beroukhimlab/Shu/ccle/poster_figs/prop_events_wide_error.pdf', width=14, height=6)
p<-ggplot(test1, aes(x=as.factor(i), y=prop_mean, fill=class)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=prop_mean-sd, ymax=prop_mean+sd), width=.2,position=position_dodge(.9))
p<- p+ scale_fill_manual(name='SV Class', values=c(alpha("#404080", 0.64), alpha("#69b3a2",0.64))) + 
   xlab('Distance to Germline Cutoff') + ylab('Proportion of Events Removed') + theme_ss
p
dev.off()





test_set<-df1

test_set[,sv_class:=ifelse(CLASS=='GERMLINE',0,1)]
test_set[,sv_pred:=ifelse(germline_dist<1000,0,1)]

test_set<-na.omit(test_set)


#graph roc/auc
pdf('/Volumes/xchip_beroukhimlab/Shu/ccle/poster_figs/roc_justgermline.pdf', width=10, height=6)
pr<-prediction(as.numeric(test_set$sv_pred), as.numeric(test_set$sv_class))
auc_ROCR <- performance(pr, measure = "auc")
auc_val <- auc_ROCR@y.values[[1]]
pref <- performance(pr, "tpr", "fpr")
plot(pref, main = paste0("ROC Curve, Germline Reference Distance Only", '\n', 'AUC: ', substr(auc_val, 1,4) ), colorize = F)
abline(a = 0, b = 1)
dev.off()




####get performance for different cutoffs 
#50,250,1000,5000

test_set_all<-df1


get_performance<-function(i, test_set){
  test_set[,sv_class:=ifelse(CLASS=='GERMLINE',0,1)]
  test_set[,sv_pred:=ifelse(germline_dist<i,0,1)]
  test_set<-na.omit(test_set)
  
  pr<-prediction(as.numeric(test_set$sv_pred), as.numeric(test_set$sv_class))
  tpr_fpr <- performance(pr, measure='tpr','fpr')
  tpr<-tpr_fpr@y.values[[1]][2]
  fpr<-tpr_fpr@x.values[[1]][2]

  
  return(data.table(tpr, fpr))
  
  
}



test1<-lapply(c(seq(0,100, 25),seq (200,900,100), seq(1000,10000, 1000)), get_performance, test_set_all)

test2<-list(1,1)
test2<-rbind(test2, rbindlist(test1))
test2<-rbind(test2, list(0, 0))

require(pracma)
AUC = trapz(test2$fpr,test2$tpr)

pdf('/Users/shu/Downloads/distance_Roc_10k.pdf')
plot(test2$fpr, test2$tpr, main = paste0("ROC Curve from Germline Reference Distance (0-10000bp)", '\n', 'AUC: ', substr(-AUC, 1,5) ), 
     xlim=c(0,1), ylim=c(0,1), xlab='False Positive Rate', ylab='True Positive Rate', colorize = F)
lines(test2$fpr, test2$tpr)
abline(a = 0, b = 1)

dev.off()



















