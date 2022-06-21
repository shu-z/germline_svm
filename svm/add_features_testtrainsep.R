library(data.table)
library(e1071)
library(caTools)
library(ROCR)
library(stringr)
library(ggfortify)
library(caret)


####AS OF 20220419
#this one to use with new annotations!!!! some columns have changed 

#read in desired dataset
#df1<-fread('/Volumes/xchip_beroukhimlab/Shu/ccle/20220425_snowman_evenlysampledtumorandsv_newannot.csv')
#this one is with evenly sampled sv as well (ish)
#df<-fread('/Volumes/xchip_beroukhimlab/Shu/ccle/20220425_snowman_evenlysampledtumor_newannot.csv')


#these datasets from 20220601 sampling (10 test tumors, 15 train tumors)
df_test<-fread('/Volumes/xchip_beroukhimlab/Shu/svm/data/20220601_snowman_newannot_test_10tumors_allsvs.csv')
df_train<-fread('/Volumes/xchip_beroukhimlab/Shu/svm/data/20220601_snowman_newannot_train_15tumors_600svs.csv')

#combine test and train
df_test[,train_test:='test']
df_train[,train_test:='train']
df<-rbind(df_test, df_train)


#encoding the target feature as factor
df$sv_class = factor(df$CLASS, levels=c('GERMLINE','SOMATIC'), labels=c(0,1))

#add a column for name (all name columns are different for germline/somatic)
split_name<-strsplit(df$name, '[.]')
sample_name<-sapply(1:nrow(df), function(i){
  return(split_name[[i]][1])
})
df$sample_name<-sample_name

#turn vcf outputs into readable features for svm
df[,homlen:= nchar(HOMSEQ)]
df[,insertion_len:= nchar(INSERTION)]
df[,hom_gc:= ifelse(nchar(HOMSEQ)>0, (str_count(HOMSEQ, 'G|C'))/nchar(HOMSEQ), 0)]
df[,insertion_gc:= ifelse(nchar(INSERTION)>0, (str_count(INSERTION, 'G|C'))/nchar(INSERTION), 0)]

#add in all svtypes and make each binary
df[,svtype := ifelse(chrom1 == chrom2, ifelse(strand1 == strand2, 'INV', ifelse(strand2 == "+", "DUP", "DEL")),"INTER")]
df[, del:=ifelse(svtype=='DEL', 1, 0)]
df[, dup:=ifelse(svtype=='DUP', 1, 0)]
df[, inv:=ifelse(svtype=='INV', 1, 0)]
df[, inter:=ifelse(svtype=='INTER', 1, 0)]

#to replace empty line/sine distances -- could be better way to do this
df[is.na(df)] <- -1

# if want to count svs/sample within df 
df[, num_sv_sample :=n_long_events]


#df_test<-df[df$train_test=='test']
#df_test<-df[df$train_test=='test']


#if we want an even somatic/germline split 
somatic_subset<-df[df$CLASS=='SOMATIC']
germline_subset<-df[df$CLASS=='GERMLINE']

#split into test and train sets 
set.seed(123)
# train_ind <- sample(seq_len(nrow(df)), size = floor(0.8 * nrow(df)))
# train <- df[train_ind, ]; test <- df[-train_ind, ]
# 


# feature scaling in train and test sets
features_toscale<-c('homlen', 'insertion_len', 'SPAN', 'gnomad_dist', 'del', 'dup', 'inv', 'inter',
                    'hom_gc', 'insertion_gc', 'line_dist', 'sine_dist', 'num_sv_sample', 'CN_annot', 'exon_annot')
train_scaled<-(train[, lapply(.SD, scale), .SDcols = features_toscale])
train_scaled<-cbind(train_scaled, sv_class=train$sv_class)
test_scaled<-(test[, lapply(.SD, scale), .SDcols = features_toscale])
test_scaled<-cbind(test_scaled, sv_class=test$sv_class)

