library(data.table)
library(stringr)



tsv_dir<-'/Users/shu/germline_svm/data/terra/'
tsvs<-list.files(tsv_dir)
tsv_pair<-tsvs[grepl('pair', tsvs)]
tsv_sample<-tsvs[grepl('sample', tsvs)]


combine_samples<-function(file_name){
  split_name<-str_split(file_name, '_')
  tt<-split_name[[1]][2]
  df<-fread(paste0(tsv_dir, file_name))
  if ('WGS_bai_path' %in% colnames(df)){
    df<-df[,c('entity:sample_id', 'participant', 'tcga_sample_id', 'sample_type',
              'WGS_bai_path', 'WGS_bam_analysis_id', 'WGS_bam_path' )]
    df<-df[!(is.na(df$WGS_bam_path) | df$WGS_bam_path==""), ]
    return(cbind(df, tumor_type=substr(tt,1, nchar(tt)-4)))
  }
}

combine_pairs<-function(file_name){
  split_name<-str_split(file_name, '_')
  tt<-split_name[[1]][2]
  df<-fread(paste0(tsv_dir, file_name))
  return(cbind(df, tumor_type=substr(tt,1, nchar(tt)-4)))
}

sample_df<-rbindlist(lapply(tsv_sample, combine_samples))
pair_df<-rbindlist(lapply(tsv_pair, combine_pairs))
pair_df<-pair_df[pair_df$participant %in% sample_df$participant]


match_pair_sample<-function(i){
  row<-pair_df[i]
  
  pair_name<-row$`entity:pair_id`
  splt_pair_name<-strsplit(pair_name, '-')

  sample_rows<-sample_df[sample_df$participant %in% row$participant]
  #print(sample_rows)
  control_row<-sample_rows[grepl(splt_pair_name[[1]][5], sample_rows$sample_type)]
  colnames(control_row) <- paste("normal", colnames(control_row), sep = "_")
  case_row<-sample_rows[grepl(splt_pair_name[[1]][4], sample_rows$sample_type)]
  colnames(case_row) <- paste("case", colnames(case_row), sep = "_")
  
  
  return(cbind(row, control_row[,4:7], case_row[,4:7]))
  
}


final_tsv<-rbindlist(lapply(1:nrow(pair_df), match_pair_sample))
final_tsv<-final_tsv[!(is.na(final_tsv$case_WGS_bam_path) | final_tsv$case_WGS_bam_path==""), ]
final_tsv<-final_tsv[!(is.na(final_tsv$normal_WGS_bam_path) | final_tsv$normal_WGS_bam_path==""), ]




#filter for equal tumor type representation
even_tt_sampling<-function(tt){
  df_tt<-final_tsv[final_tsv$tumor_type == tt]
  return(df_tt[sample(seq_len(nrow(df_tt)), size = 5),])
}

sampled_tsv<-rbindlist(lapply(unique(final_tsv$tumor_type), even_tt_sampling))
fwrite(sampled_tsv, paste0(tsv_dir, 'tcga_sampled_wgs.tsv'), sep='\t')


#do only dlbcl for now
dlbcl_tsv<-final_tsv[final_tsv$tumor_type=='dlbc']
fwrite(dlbcl_tsv, paste0(tsv_dir, 'dlbcl_wgs.tsv'), sep='\t')


