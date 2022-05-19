
#############################################
# additional graphs



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
test_df<-cbind(test, pred=(y_pred), prob=max_prob)
test_df[,missclassified:=ifelse(sv_class==y_pred, 0, 1)]

#take substr of sample name 
split_name<-strsplit(test_df$sample, '[.]')
short_name<-lapply(1:length(split_name), function(i){
  return(split_name[[i]][1])
})
test_df$patient_id<-short_name


test_df_germline<-test_df[test_df$CLASS=='GERMLINE']

missclass_by_sample<-function(sample_name){
  samp<-test_df[test_df$patient_id==sample_name]
  class_correct<-length(which(samp$missclassified==0))
  return(data.table(prop_class=(class_correct/nrow(samp)),
          tumor_type=substr(samp$project_code[1], 1, nchar(samp$project_code[1])-3)))
}

count_missclass<-rbindlist(lapply(unique(test_df$patient_id), missclass_by_sample))
count_missclass_germline<-rbindlist(lapply(unique(test_df_germline$patient_id), missclass_by_sample))


p<-ggplot(count_missclass, aes(x=tumor_type, y=prop_class, fill=tumor_type)) +
  geom_boxplot() + geom_point(aes(), size = 1, shape = 21) + ylim(c(0,1)) + 
  theme_ss
# p<-p + labs(title='Proportion of Germline SVs Removed', subtitle='16k train, 4k test',
#             x='tumor type', y='proportion of events removed')
# p<-p + labs(title='Proportion of Correctly Classified Events', subtitle='16k train, 4k test',
#  x='Tumor Type', y='Proportion of Events Correctly Classified')

p<-p + labs(y='Proportion of Events Correctly Classified', x='')  +  
    guides(fill=guide_legend(title="Tumor Type"))
pdf('/Users/shu/germline_svm/figs/20220321_boxplot_events_classified.pdf', width=10, height=5)
p
dev.off()






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