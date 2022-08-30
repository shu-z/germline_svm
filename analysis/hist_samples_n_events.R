library(data.table)
library(stringr)
library(ggfortify)
library(ggplot2)



#SCRIPT TO LOOK AT DISTRIBUTION OF SVS WITH NEW ANOTATION SITUATION 

pcawg_metadata<-fread('/Volumes/xchip_beroukhimlab/Alex/PCAWG/PCAWG_metadata_OCT2020.csv', sep=',')


#function to map sample name to type 
map_id<-function(i){
  sample_df<-fread(i, nrows=1)
  sample_name<-substring(sample_df$sample, 1, 36)
  project_code<-pcawg_metadata$dcc_project_code[pcawg_metadata$tumor_wgs_submitter_sample_id %in% sample_name]
  return(data.table(cbind(path=i, sample_name, project_code)))
}

#read in snowman germline somatic calls 
#this is updated one as of 04/19/2022, with updated coding annotations
snowman_germdist_dir<-'/Volumes/xchip_beroukhimlab/Alex/CCLE_PCAWG/Figures_code/figure_data/data/4_annot/'

all_paths<-paste0(snowman_germdist_dir, list.files(snowman_germdist_dir))
map_df<-rbindlist(lapply(all_paths, map_id), fil=T)
#sort(table(map_df$project_code))   #min num tumors is 7 from dlbcl 


#remove short events and select same number of svs 
even_sv_sampling<-function(row){
  sample_df<-fread(row[1])
  long_sample_df<-sample_df[(sample_df$SPAN>=1e3) | (sample_df$SPAN ==-1)]
  n_somatic<-sum(long_sample_df$CLASS=='SOMATIC')
  n_germline<-sum(long_sample_df$CLASS=='GERMLINE')
  return(data.table(cbind(sample=row[2], project_code=row[3], n_long=nrow(long_sample_df),  n_somatic, n_germline)))
}

final_svm_df<-rbindlist(apply(map_df, 1, even_sv_sampling))
#write.csv(final_svm_df, '/Users/shu/germline_svm/data/20220526_updatedannotations_longevents_persample.csv')
final_svm_df<-fread('/Users/shu/germline_svm/data/20220526_updatedannotations_longevents_persample.csv', drop=1)

final_svm_df[, n_long:=as.numeric(n_long)]
final_svm_df[, n_somatic:=as.numeric(n_somatic)]
final_svm_df[, n_germline:=as.numeric(n_germline)]
final_svm_df[, prop_somatic:=n_somatic/n_long]


mean_propsomatic<-aggregate(prop_somatic ~ project_code, final_svm_df, 'mean')


#if want bins in log space 
log_space<-unique(round(exp(seq(log(1), log(35000), length.out=50))))



hist(as.numeric(final_svm_df$n_long), 1200, xlim=c(0,3000),
     xlab='# Long SVs per Sample', ylab='Count', main='')



#ggplot theme 
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



#make boxplot to look at #SVs per tumor type, and soomatic/germline proportion 

p<-ggplot(final_svm_df, aes(x=project_code, y=n_long, fill=project_code)) +
  geom_boxplot() + geom_point(aes(), size = 1, shape = 21) + 
  scale_y_continuous(trans='log10') + 
  #geom_violin() + geom_point(aes(), size = 1, shape = 21) + ylim(c(0,1)) + 
  theme_ss
p<-p + labs(title='Long SVs per Sample (SVM training)', x='Tumor Type', y='# Long SVs per Sample')  +  
  guides(fill=guide_legend(title="Tumor Type"))
pdf('/Users/shu/germline_svm/figs/20220526_boxplot_longevents_count.pdf', width=10, height=5)
p
dev.off()


