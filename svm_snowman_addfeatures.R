library(data.table)
library(e1071)
library(caTools)
library(ROCR)
library(stringr)
library(ggfortify)
library(caret)


#########################################
# code to filter events >1k (don't need to run every time)
# aims to select equal # of SVs from equal # of tumor types 

pcawg_metadata<-fread('/Volumes/xchip_beroukhimlab/Alex/PCAWG/PCAWG_metadata_OCT2020.csv', sep=',')
#pcawg_metadata<-fread('/Users/shu/ClusterSV/supp_table_pcawg.csv', sep=',')
#which(pcawg_metadata == '00493087-9d9d-40ca-86d5-936f1b951c93', arr.ind=T)

#function to map sample name to type 
map_id<-function(i){
  sample_df<-fread(i, nrows=1)
  sample_name<-substring(sample_df$sample, 1, 36)
  project_code<-pcawg_metadata$dcc_project_code[pcawg_metadata$tumor_wgs_submitter_sample_id %in% sample_name]
  return(data.table(cbind(path=i, sample_name, project_code)))
}

#read in snowman germline somatic calls 
#new one has line-sine annotation
snowman_germdist_dir<-'/Volumes/xchip_beroukhimlab/Alex/CCLE_PCAWG/Germline_Correctiion/Test_Set/Germ_dist/FUNCT_ANNOT/Line_Sine_Annot/'
all_paths<-paste0(snowman_germdist_dir, list.files(snowman_germdist_dir))
map_df<-rbindlist(lapply(all_paths, map_id), fil=T)
#table(map_df$project_code)   #min num tumors is 7 from dlbcl 

#select sample number of files per tutmor type 
even_tt_sampling<-function(tumor_type){
  tumor_paths<-map_df$path[map_df$project_code == tumor_type]
  return(data.table(paths=sample(tumor_paths, size=7), project_code=tumor_type))
}

#remove short events and select same number of svs 
even_sv_sampling<-function(row){
  sample_df<-fread(row[1])
  long_sample_df<-sample_df[(sample_df$SPAN>=1e3) | (sample_df$SPAN ==-1)]
  print(nrow(long_sample_df))
  n_len<-min(nrow(long_sample_df), 1200)
  return(cbind(long_sample_df[sample(seq_len(nrow(long_sample_df)), size = n_len),], 
               n_long_events=nrow(long_sample_df), project_code=row[2], file_path=row[1]))
}
  
sample_paths<-rbindlist(lapply(unique(map_df$project_code), even_tt_sampling))
final_svm_df<-rbindlist(apply(sample_paths, 1, even_sv_sampling))
write.csv(final_svm_df, '/Volumes/xchip_beroukhimlab/Shu/ccle/20220302_snowman_evenlysampledtumor.csv')


#########################################

##run everything below here 
#read in desired dataset
#this one is filtered for SVs >1kb, equal number of tumors and SVs/tumor 
df_all<-fread('/Volumes/xchip_beroukhimlab/Shu/ccle/20220302_snowman_evenlysampled.csv')
df<-df_all

# Encoding the target feature as factor
df$sv_class = factor(df$CLASS, levels=c('GERMLINE','SOMATIC'), labels=c(0,1))

#add in features to test 
df[,homlen:= nchar(HOMSEQ)]
df[,insertion_len:= nchar(INSERTION)]
df[name == df$name, svtype := ifelse(chrom1 == chrom2, ifelse(strand1 == strand2, 
                              ifelse(strand1 == '+', "h2hINV", "t2tINV"), ifelse(strand2 == "+", "DUP", "DEL")),"INTER")]
df$sv_type_factor = as.integer(factor(df$svtype, levels=c('DUP', 'DEL', 't2tINV', 'h2hINV', 'INTER'), labels=1:5))
df[,hom_gc:= ifelse(nchar(HOMSEQ)>0, (str_count(HOMSEQ, 'G|C'))/nchar(HOMSEQ), 0)]
df[,insertion_gc:= ifelse(nchar(INSERTION)>0, (str_count(INSERTION, 'G|C'))/nchar(INSERTION), 0)]
df[,annot_class := as.integer(factor(annot, levels=c('NON_PROTEIN_CODING','INTRON', 'CODING', 'INTERCHROMOSOMAL', 'CN'), labels=c(0:4)))]
#to replace empty line/sine distances -- could be better way to do this
df[is.na(df)] <- -1

# if want to count svs/sample within df 
# sample_counts<-data.table(table(df$file_path))
# file_paths<-df$file_path
# sample_svs<-sapply(1:nrow(df), function(i){
#   return(sample_counts$N[sample_counts$V1==file_paths[i]])
# })
df[, num_sv_sample :=n_long_events]


#if we want an even somatic/germline split 
somatic_subset<-df[df$CLASS=='SOMATIC']
germline_subset<-df[df$CLASS=='GERMLINE']
df<-rbind(somatic_subset[sample(seq_len(nrow(somatic_subset)), size = min(nrow(somatic_subset),10000)),],
          germline_subset[sample(seq_len(nrow(germline_subset)), size = min(nrow(germline_subset), 10000)),])

#split into test and train sets 
set.seed(123)
train_ind <- sample(seq_len(nrow(df)), size = floor(0.8 * nrow(df)))
train <- df[train_ind, ]; test <- df[-train_ind, ]


# feature scaling
features_toscale<-c('homlen', 'insertion_len', 'SPAN', 'germline_dist', 'sv_type_factor', 
                    'hom_gc', 'insertion_gc', 'line_dist', 'sine_dist', 'annot_class')
train_scaled<-(train[, lapply(.SD, scale), .SDcols = features_toscale])
train_scaled<-cbind(train_scaled, sv_class=train$sv_class)
test_scaled<-(test[, lapply(.SD, scale), .SDcols = features_toscale])
test_scaled<-cbind(test_scaled, sv_class=test$sv_class)



#use radial if not linear kernel
classifier = svm(formula = sv_class ~ .,
                 data = train_scaled, type = 'C-classification',
                 kernel = 'radial', cachesize=400, 
                 probability=T)

#predict test set with classifier 
y_pred <- predict(classifier, newdata = test_scaled, decision.values = T, probability = T)
probabilities<-attr(y_pred, 'probabilities')

max_prob<-apply(probabilities,1, function(x){
  val<-as.vector(x)
  #print(val)
  test<-which.max(val)-1
  ifelse(test==0, return(val[2]), return(val[2]))
})

test_scaled<-na.omit(test_scaled)


#graph roc/auc
#pdf('/Volumes/xchip_beroukhimlab/Shu/ccle/poster_figs/20211214_roc_50k_10features_wide.pdf', width=10, height=6)
#pr<-prediction(as.numeric(y_pred), as.numeric(test_scaled$sv_class))
pr<-prediction(as.numeric(max_prob), as.numeric(test_scaled$sv_class))

auc_ROCR <- performance(pr, measure = "auc")
auc_val <- auc_ROCR@y.values[[1]]
pref <- performance(pr, "tpr", "fpr")
#plot(pref, main = paste0("ROC curve", '\n', 'AUC: ', auc_val), colorize = F)
plot(pref, main = paste0('AUC: ', substring(auc_val,1,5) ), colorize = F, 
     cex.lab=1.5, cex.main=1.5)

abline(a = 0, b = 1)
#dev.off()


#tune classifier 
obj <- tune(svm, sv_class~., data = train_scaled, 
            #ranges = list(gamma = 2^(-6:2), cost = 2^(0:4)),
            ranges = list(gamma = 2^(-1:2), cost = 10^(1:5)),
           tunecontrol = tune.control(sampling = "fix"))
summary(obj)

saveRDS(obj, '/Volumes/xchip_beroukhimlab/Shu/ccle/poster_figs/svm_tuning_obj.rds')
pdf('/Volumes/xchip_beroukhimlab/Shu/ccle/poster_figs/svm_performance.pdf')
plot(obj)
dev.off()


##### try platt's calibration
source('/Users/shu/GermlineSVAnnotator/svm/prcalibrate_eric.R')

test1<-prCalibrate(as.numeric(levels(test_scaled$sv_class))[test_scaled$sv_class], as.numeric(max_prob))
pr_cal<-prediction((test1$cal.probs), as.numeric(test_scaled$sv_class))

auc_ROCR_cal <- performance(pr_cal, measure = "auc")
auc_val_cal <- auc_ROCR_cal@y.values[[1]]
pref_cal <- performance(pr_cal, "tpr", "fpr")
#plot(pref, main = paste0("ROC curve", '\n', 'AUC: ', auc_val), colorize = F)
plot(pref_cal, main = paste0('AUC: ', substring(auc_val_cal,1,5) ), colorize = F, 
     cex.lab=1.5, cex.main=1.5)

abline(a = 0, b = 1)





#############################################
# additional graphs

###### look at pca of variables 
pca_df<-na.omit(test_scaled)
pca_df <- pca_df[is.finite(rowSums(pca_df)),]
pc <- prcomp(pca_df[,1:10],
             center = TRUE,
             scale. = TRUE)
attributes(pc)

#pdf('/Volumes/xchip_beroukhimlab/Shu/ccle/poster_figs/20211209_pca_nolog.pdf')
p<-autoplot(pc, data = pca_df, colour = 'sv_class',  
            loadings = TRUE, loadings.colour = 'blue',
            loadings.label = TRUE, loadings.label.size = 3,loadings.label.repel=T)
p<-p+ geom_point(alpha = 0.01)+ theme_minimal() 
p
dev.off()

train_scaled<-na.omit(train)
train_scaled$sv_class<-as.numeric(train_scaled$sv_class)

### plotting boundaries for every pair of features 
pdf('/Volumes/xchip_beroukhimlab/Shu/ccle/20211213_svm_classify_allvar_smaller.pdf', width=4, height=4)
all_indices<-data.table(expand.grid(1:10, 1:10))
plot_func<-function(i, train_scaled, classifier){
  a<-all_indices[[i,1]]
  b<-all_indices[[i,2]]
  if (a>b){
    test<-colnames(train_scaled)
    plot(classifier, train_scaled, as.formula(paste(test[[b]], '~', test[[a]])),
         fill = TRUE, grid=10, symbolPalette = c("#E30707",'#0DEE58'),
         #xlim=c(-1,10), ylim=c(-1,10),
         col=c( alpha("#404080", 0.6), alpha("#69b3a2", 0.6))
    )
  }
}

lapply(1:nrow(all_indices), plot_func, train_scaled, classifier)
dev.off()

### if want to plot a single pair of features
#pdf('/Volumes/xchip_beroukhimlab/Shu/ccle/20211118_svm_classify_germlinedist_span_zoom1.pdf')
plot(classifier, train_scaled, germline_dist.V1 ~ SPAN.V1, 
     fill = TRUE, grid=10, symbolPalette = c("#E30707",'#0DEE58'),
     xlim=c(-1,3), ylim=c(-1,3),
     #color.palette = terrain.colors
     
     col=c( alpha("#404080", 0.6), alpha("#69b3a2", 0.6))
     #col=cols
)
#dev.off()


### plot histogram features
scale<-2
theme_ss <- theme_classic() +
  #theme(axis.text.x = element_text(angle = 45, vjust = -0.5, size=8*scale),
  theme(axis.text.x = element_text( size=16*scale),
        axis.text.y = element_text( size=12*scale),
        axis.title=element_text(size=21*scale),
        plot.margin =  unit(c(0.5,0,1.5,0), "cm") )



features_hist<-c(  'germline_dist', 'line_dist', 'sine_dist', 'homlen', 'hom_gc', 'SPAN',
                   'insertion_len', 'svtype', 'annot', 'sv_class')
# features_hist<-c(  'germline_dist', 'line_dist', 'sine_dist', 'homlen', 'hom_gc', 'SPAN',
#                     'svtype', 'annot', 'sv_class')


to_plot<-train[,..features_hist]
to_plot[,annot:=gsub('NON_PROTEIN_CODING', 'NONCODING', annot)]
to_plot[,annot:=gsub('INTERCHROMOSOMAL', 'INTER', annot)]
to_plot[,annot:=gsub("^$", NA, annot)]

#to_plot$annot = with(to_plot, reorder(class, hwy, median))



new_name<-c( 'Distance to Reference Germline SV', 'Distance to LINEs', 'Distance to SINEs',
             'Breakpoint Homology Length','Homology GC','SV Span',
             'Insertion Length', 'Type of SV', 'Breakpoint Annotation')

plot_hist<-function(i){
  
  p <- to_plot %>%
    ggplot( aes_string(x=(features_hist[i]), fill='sv_class')) +
    geom_histogram( color="#e9ecef", alpha=0.7, position = 'identity', bins=100) +
    scale_fill_manual(values=c( "#404080", "#69b3a2")) +
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))+
    theme_ss + labs(fill="") + xlab(new_name[i]) + ylab('Count') + theme(legend.position = "none")
}
plots1<-lapply(1:6, plot_hist)

plot_bar<-function(i){
  
  p<-ggplot(subset(to_plot, !is.na(annot)), aes_string(x=features_hist[i], fill='sv_class')) +
    geom_bar(position='dodge', alpha=0.6) + xlab(new_name[i]) +ylab('Count') + 
    scale_fill_manual(values=c( "#404080", "#69b3a2")) +theme_ss + theme(legend.position = "none")
  p
}
plots2<-lapply(8:9, plot_bar)

pdf('/Volumes/xchip_beroukhimlab/Shu/ccle/poster_figs/20211214_featurehist_8features_5520.pdf', width=55, height=25)
cowplot::plot_grid(plotlist = c(plots1[1:3],plots2[1], plots1[4:6], plots2[2]), nrow = 2)

dev.off()



#### look at events that were misclassified 
events_mis<-cbind(na.omit(test), y_pred)
events_mis<-events_mis[events_mis$sv_class!=events_mis$y_pred]

events_mis[,misclass:=ifelse(sv_class==0, 'germline_misclassified_somatic', 'somatic_misclassfied_germline')]

#germline_class_somatic<-events_mis[events_mis$sv_class==0]
#somatic_class_germline<-events_mis[events_mis$sv_class==1]


misclass_hist<-c('germline_dist', 'line_dist', 'sine_dist', 'homlen', 'hom_gc', 'SPAN',
                 'insertion_len', 'svtype', 'annot', 'sv_class', 'misclass')


plot_misclass<-events_mis[,..misclass_hist]
plot_misclass[,annot:=gsub('NON_PROTEIN_CODING', 'NONCODING', annot)]
plot_misclass[,annot:=gsub('INTERCHROMOSOMAL', 'INTER', annot)]
plot_misclass[,annot:=gsub("^$", NA, annot)]


new_name<-c( 'Distance to Reference Germline SV', 'Distance to LINEs', 'Distance to SINEs',
             'Breakpoint Homology Length','Homology GC','SV Span',
             'Insertion Length', 'Type of SV', 'Breakpoint Annotation')

plot_hist_misclass<-function(i){
  
  p <- plot_misclass %>%
    ggplot( aes_string(x=(misclass_hist[i]), fill='misclass')) +
    geom_histogram( color="#e9ecef", alpha=0.7, position = 'identity', bins=150) +
    scale_fill_manual(values=c( "#88B159", "#df6b39")) +
    scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))+
    theme_ss + labs(fill="") + xlab(new_name[i]) + ylab('Count') + theme(legend.position = "none")
  print(p)
}
plots1<-lapply(1, plot_hist_misclass)

plot_bar_misclass<-function(i){
  
  p<-ggplot(subset(plot_misclass, !is.na(annot)), aes_string(x=features_hist[i], fill='misclass')) +
    geom_bar(position = position_dodge(preserve = "single"), alpha=0.5) + xlab(new_name[i]) +ylab('Count') + 
    scale_fill_manual(values=c( "#88B159", "#df6b39")) +theme_ss + theme(legend.position = "none")
  print(p)
}
plots2<-lapply(8:9, plot_bar_misclass)

pdf('/Volumes/xchip_beroukhimlab/Shu/ccle/poster_figs/20211214_features_missclass.pdf', width=20, height=18)
#cowplot::plot_grid(plotlist = c(plots1[1:4],plots2[1], plots1[6], plots2[2]), nrow = 2)
cowplot::plot_grid(plotlist = c(plots1[1],plots2[1]), nrow = 2)

dev.off()

























