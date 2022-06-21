library(data.table)

#############################################
# additional graphs to look at misclassified events 



scale<-1.5
theme_ss <- theme_bw(base_size=16) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.text.x = element_text(angle = 90, vjust = 0.5, size=8*scale, family="mono"),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y = element_text(hjust = 0.5,size=10*scale)
        #axis.text = element_text(size = 10*scale, family = "mono"))
  )



#make boxplot for correctly classified events by tissue type 
test_misclass<-cbind(test, pred=(y_pred), prob=as.numeric(probabilities$`1`))

test_misclass[,missclassified:=ifelse(sv_class==y_pred, 0, 1)]

#take substr of sample name 
split_name<-strsplit(test_misclass$sample, '[.]')
short_name<-lapply(1:length(split_name), function(i){
  return(split_name[[i]][1])
})
test_misclass$patient_id<-short_name


test_misclass_germline<-test_misclass[test_misclass$CLASS=='GERMLINE']
test_misclass_somatic<-test_misclass[test_misclass$CLASS=='SOMATIC']

ggplot(germline_missclass, aes(gnomad_dist)) +              
  geom_histogram(bins = 100) +
  scale_x_log10() + theme_classic()

missclass_by_sample<-function(pat_id){
  samp<-test_misclass[test_misclass$patient_id==pat_id]
  class_correct<-length(which(samp$missclassified==0))
  return(data.table(prop_class=(class_correct/nrow(samp)),
                    tumor_type=substr(samp$project_code[1], 1, nchar(samp$project_code[1])-3)))
}

count_missclass<-rbindlist(lapply(unique(test_misclass$patient_id), missclass_by_sample))
count_missclass_germline<-rbindlist(lapply(unique(test_misclass_germline$patient_id), missclass_by_sample))


p<-ggplot(count_missclass, aes(x=tumor_type, y=prop_class, fill=tumor_type)) +
  geom_boxplot() + geom_point(aes(), size = 1, shape = 21) + ylim(c(0,1)) + 
  theme_ss
# p<-p + labs(title='Proportion of Germline SVs Removed', subtitle='16k train, 4k test',
#             x='tumor type', y='proportion of events removed')
# p<-p + labs(title='Proportion of Correctly Classified Events', subtitle='16k train, 4k test',
#  x='Tumor Type', y='Proportion of Events Correctly Classified')

p<-p + labs(y='Proportion of Events Correctly Classified', x='')  +
  guides(fill=guide_legend(title="Tumor Type"))
#pdf('/Users/shu/germline_svm/figs/20220321_boxplot_events_classified.pdf', width=10, height=5)
p
dev.off()





#### look at events that were misclassified 
events_mis<-test_misclass[test_misclass$sv_class!=test_misclass$pred]

events_mis[,misclass:=ifelse(sv_class==0, 'germline_misclassified_somatic', 'somatic_misclassfied_germline')]

#germline_class_somatic<-events_mis[events_mis$sv_class==0]
#somatic_class_germline<-events_mis[events_mis$sv_class==1]


misclass_hist<-c('gnomad_dist', 'line_dist', 'sine_dist',  'SPAN', 'homlen', 'insertion_len','num_sv_sample', 
                 'del', 'dup', 'inv', 'inter', 'CN_annot', 'exon_annot', 'sv_class', 'misclass')

plot_misclass<-events_mis[,..misclass_hist]



new_name<-c( 'Distance to Reference Germline SV', 'Distance to LINEs', 'Distance to SINEs','SV Span',
             'Breakpoint Homology Length','Insertion Length', 'Num Long SVs per Sample', 
             'Deletion', 'Duplication', 'Inversion', 'Interchromosomal', 'CN Annot', 'Exon Annot')

plot_hist_misclass<-function(i){
  
  p <- plot_misclass %>%
    ggplot( aes_string(x=(misclass_hist[i]), fill='misclass')) +
    geom_histogram( color="#e9ecef", alpha=0.85, position = 'identity', bins=150) +
    scale_fill_manual(values=c( "#88B159", "#df6b39")) +
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))+
    theme_ss + labs(fill="") + xlab(new_name[i]) + ylab('Count') + theme(legend.position = "none")
  print(p)
}
#pdf('/Users/shu/germline_svm/figs/20220615_missclassified_continuousvars.pdf', width=8, height=5)
plots1<-lapply(1:7, plot_hist_misclass)
#dev.off()




#sv class df
svtype_var<-c('svtype', 'sv_class', 'pred')
svtype_df<-test_misclass[, ..svtype_var]
svtype_df[,missclassified:=ifelse(sv_class==y_pred, 0, 1)]
svtype_df[,misclass:=ifelse(missclassified==0, ifelse(sv_class==0, 'Correct Germline Classification', 'Correct Somatic Classification'),
                      ifelse(sv_class==0, 'Misclassified Germline', 'Misclassified Somatic'))]

# p <- ggplot( svtype_df, aes(x=svtype, fill=misclass)) +
#   geom_bar(position = position_dodge(), alpha=1) +ylab('Count') + 
#   scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))+
#   xlab('SV Type') + 
#   scale_fill_manual(values=c('#F8766D', '#00BFC4', "#88B159", "#df6b39")) +theme_ss 
# pdf('/Users/shu/germline_svm/figs/20220615_missclassified_svtype_4bar.pdf', width=10, height=6)
# (p)
# dev.off()


#count up occurrences in r 
count_occurrences<-function(svt){
  svtype1_df<-svtype_df[svtype_df$svtype==svt]
  prop_germ_misclass<-sum(svtype1_df$misclass=='Misclassified Germline')/sum(svtype1_df$sv_class==0)
  prop_som_misclass<-sum(svtype1_df$misclass=='Misclassified Somatic')/sum(svtype1_df$sv_class==1)
  
  return(data.table('SV_Type'=svt,Proportion= c(prop_germ_misclass, prop_som_misclass), 
                    Misclassification=c('Germline Misclass', 'Somatic Misclass')))
}

count_svtypes<-rbindlist(lapply(unique(svtype_df$svtype), count_occurrences))



p<-ggplot(count_svtypes, aes(SV_Type, Proportion)) +   
  geom_bar(aes(fill = Misclassification), position = "dodge", stat="identity") + 
  xlab('SV Type') + ylab ('Proportion of Misclassified SVs') + theme_ss + 
  theme(axis.text.x = element_text(vjust = 0.5, size=8*scale)) + 
  scale_fill_manual(values=c( "#88B159", "#df6b39"))  


pdf('/Users/shu/germline_svm/figs/20220615_missclassified_svtype.pdf', width=10, height=6)
(p)
dev.off()



annot_var<-c('CN_annot', 'exon_annot', 'sv_class', 'pred')
annot_df<-test_misclass[, ..annot_var]
annot_df[,missclassified:=ifelse(sv_class==y_pred, 0, 1)]
annot_df[,misclass:=ifelse(missclassified==0, ifelse(sv_class==0, 'Correct Germline Classification', 'Correct Somatic Classification'),
                            ifelse(sv_class==0, 'Misclassified Germline', 'Misclassified Somatic'))]

count_occurrences_annot<-function(annot){
  annot_col<-c(annot, 'sv_class', 'pred', 'missclassified', 'misclass')
  annot1_df<-annot_df[, ..annot_col]
  print(annot1_df[1:5])
  prop_germ_misclass<-sum(annot1_df$misclass=='Misclassified Germline')/sum(annot1_df[,1]==0)
  prop_som_misclass<-sum(annot1_df$misclass=='Misclassified Somatic')/sum(annot1_df[,1]==1)
  
  return(data.table('Annot'=annot,Proportion= c(prop_germ_misclass, prop_som_misclass), 
                    Misclassification=c('Germline Misclass', 'Somatic Misclass')))
}

count_annot<-rbindlist(lapply(c('CN_annot', 'exon_annot'), count_occurrences_annot))



p<-ggplot(count_annot, aes(Annot, Proportion)) +   
  geom_bar(aes(fill = Misclassification), position = "dodge", stat="identity") + 
  xlab('SV Type') + ylab ('Proportion of Misclassified SVs') + theme_ss + 
  theme(axis.text.x = element_text(vjust = 0.5, size=8*scale)) + 
  scale_fill_manual(values=c( "#88B159", "#df6b39"))  


pdf('/Users/shu/germline_svm/figs/20220615_missclassified_annot.pdf', width=10, height=6)
(p)
dev.off()



#pdf('/Volumes/xchip_beroukhimlab/Shu/ccle/poster_figs/20211214_features_missclass.pdf', width=20, height=18)
#cowplot::plot_grid(plotlist = c(plots1[1:4],plots2[1], plots1[6], plots2[2]), nrow = 2)
#cowplot::plot_grid(plotlist = c(plots1[1],plots2[1]), nrow = 2)

#dev.off()


#######cluster misclassified events 
library(factoextra)

test_misclass_scaled<-cbind(test_scaled, pred=(y_pred), prob=as.numeric(probabilities$`1`))
test_misclass_scaled[,missclassified:=ifelse(sv_class==y_pred, 0, 1)]
events_mis_scaled<-test_misclass_scaled[test_misclass_scaled$sv_class!=test_misclass_scaled$pred]
events_mis_scaled[,sv_class:=as.numeric(sv_class)]
events_mis_scaled[,pred:=as.numeric(pred)]

#events_mis[,misclass:=ifelse(sv_class==0, 'germline_misclassified_somatic', 'somatic_misclassfied_germline')]

test_feat<-cbind(test[, ..features_toscale], sv_class=test$sv_class)
test_misclass<-cbind(test_feat, pred=(y_pred), prob=as.numeric(probabilities$`1`))
test_misclass[,missclassified:=ifelse(sv_class==y_pred, 0, 1)]
events_mis<-test_misclass[test_misclass$sv_class!=test_misclass$pred]
events_mis[,sv_class:=as.numeric(sv_class)]
events_mis[,pred:=as.numeric(pred)]
#events_mis_scaled<-scale(events_mis[,1:16])

events_mis_germline<-events_mis_scaled[events_mis_scaled$sv_class==1]
events_mis_somatic<-events_mis_scaled[events_mis_scaled$sv_class==2]

fviz_nbclust(events_mis_scaled, kmeans, method = "wss")


events_mis_test<-events_mis_germline
wss <- function(k) {
  kmeans(events_mis_test, k, nstart = 10 )$tot.withinss
}

k.values<-1:20
# extract wss for 2-15 clusters

wss_values <- map_dbl(k.values, wss)
plot(k.values, wss_values, type="b", pch = 19, cex=0.4,
     xlab="Number of clusters K", ylab="Total within-clusters sum of squares")



fit <- kmeans(events_mis_test, 4) # 5 cluster solution
# get cluster means
aggregate(events_mis_test,by=list(fit$cluster),FUN=mean)
# append cluster assignment
events_mis_test <- data.frame(events_mis_test, fit$cluster)


fviz_cluster(fit, data = events_mis_test, geom = "point", scale=F)+ ggtitle("k = 4")

