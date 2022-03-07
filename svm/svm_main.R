library(data.table)
library(e1071)
library(caTools)
library(ROCR)
library(stringr)
library(ggfortify)
library(caret)


#read in desired dataset
#this one is filtered for SVs >1kb, equal number of tumors and SVs/tumor 
df<-fread('/Volumes/xchip_beroukhimlab/Shu/ccle/20220302_snowman_evenlysampled.csv')


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
  return(val[2])
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

#saveRDS(obj, '/Volumes/xchip_beroukhimlab/Shu/ccle/poster_figs/svm_tuning_obj.rds')
#pdf('/Volumes/xchip_beroukhimlab/Shu/ccle/poster_figs/svm_performance.pdf')
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










############################
#MAIN
option_list <- list(
  
  make_option(c("-i", "--input"),  type = "character", default = NULL,  help = "Input to annotated SVs"),
  make_option(c("-s", "--sort"),  type = "character", default = NULL,  help = "Indicates if SVs needs to be sorted"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "Output directory path"),
)



# inputs
parseobj = OptionParser(option_list = option_list)
opt = parse_args(parseobj)

if (!file.exists(opt$input)) {
  stop(sprintf("Input file '%s' does not exist! ", opt$input))
} 

else {
  input_pth = paste0(opt$input)
  sample <- unlist(strsplit(input_pth, "/"))[length(unlist(strsplit(input_pth, "/")))]
  #sample <- gsub(".vcf.bedpe","", sample)
  out_pth = paste0(opt$output)
  
  cat('Building references... \n')
  
  # get reference gene annotations and germline
  hg38_annotgenes = fread(paste0(opt$genomeant))
  hg38_exon_locations = fread(paste0(opt$exonannot))
  hg38_germline_gnomad <<- fread(paste0(opt$germline))
  
  
  ### build granges for reference
  gene_locations <<- GRanges(hg38_annotgenes$chromosome, IRanges(hg38_annotgenes$start, 
                                                                 hg38_annotgenes$end), EnsemblID = hg38_annotgenes$Ensembl_ID, GeneName = hg38_annotgenes$gene_name)
  
  exon_locations <<- GRanges(hg38_exon_locations$Chromosome, IRanges(hg38_exon_locations$Exon_start, 
                                                                     hg38_exon_locations$Exon_end), Exon_Ensembl_ID = hg38_exon_locations$Exon_Ensembl_ID)
  
  cat("Reading file...\n")
  bedpe_inp <- fread(paste0(input_pth))
  
  cat("Adding gene annotations...\n")
  bedpe_gene_annotation <- rbindlist(mclapply(1:nrow(bedpe_inp), annotate_sv, bedpe_inp, mc.cores = opt$cores))
  
  cat("Adding exon annotations...\n")
  bedpe_gene_exon_annotation <- rbindlist(mclapply(1:nrow(bedpe_gene_annotation), annotate_exon, bedpe_gene_annotation, mc.cores = opt$cores))
  
  cat("Fuzzy filtering germline...\n")
  bedpe_fuzzy_filtered <- rbindlist(mclapply(1:nrow(bedpe_gene_exon_annotation), fuzzy_filter_germline, bedpe_gene_exon_annotation, mc.cores = opt$cores))
  bedpe_fuzzy_filtered_sub <- bedpe_fuzzy_filtered[,uuid:=NULL]
  bedpe_somatic_only <- bedpe_fuzzy_filtered_sub[grep("Somatic", Filter)]
  bedpe_somatic_only[,Filter := NULL]
  output_somatic_only <- paste0(out_pth, sample, "_somatic_only_sv_annotated.bedpe")
  output_all <- paste0(out_pth, sample, "_sv_annotated.bedpe")
  
  write.table(bedpe_fuzzy_filtered_sub, output_all, sep = '\t', row.names = F, col.names = T, quote = F)
  write.table(bedpe_somatic_only, output_somatic_only, sep = '\t', row.names = F, col.names = T, quote = F)
  
  time_end <- Sys.time()
  cat(paste0("Began at ", time_begin,"\n"))
  cat(paste0("Ended at ", time_end,"\n"))
  
}