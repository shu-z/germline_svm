library(data.table)

#additional analyses 



germline<-fread('/Users/Shu/germline_svm/data/20220511_snowman_allsamples_newannot_germlineonly.csv')


#to do straight from reading in stuff
get_svprops<-function(path){
  df1<-fread(path)
  df<-df1[(df1$SPAN>=50) | (df1$SPAN ==-1)]
  df<-df[df$CLASS == 'SOMATIC']
  df[,svtype := ifelse(chrom1 == chrom2, ifelse(strand1 == strand2, 'INV', ifelse(strand2 == "+", "DUP", "DEL")),"INTER")]
  sample_name<-strsplit(df$name, '[.]')[[1]][1]
  
  project_code<-pcawg_metadata$dcc_project_code[pcawg_metadata$tumor_wgs_submitter_sample_id %in% sample_name]
  sv_tb<-table(df$svtype)
  return(data.table(sample_name, project_code, del=sv_tb['DEL'], dup=sv_tb['DUP'], 
                    inter=sv_tb['INTER'], inv=sv_tb['INV'], sv_sum=sum(sv_tb)))
  
  
}


sv_types_df<-lapply(all_paths, get_svprops)
sv_types_df1<-rbindlist(sv_types_df)
#read in df of all svs

#df<-fread('/Volumes/xchip_beroukhimlab/Shu/ccle/20220511_snowman_allsamples_newannot.csv')
df<-fread('/Volumes/xchip_beroukhimlab/Shu/ccle/202200419_snowman_259samples_newannot.csv')


#see if sv type proportions different across lineage
get_svtypes_pertumor<-function(i){
  df_samp<-somatic_subset[somatic_subset$sample_name==i]
  project_code<-pcawg_metadata$dcc_project_code[pcawg_metadata$tumor_wgs_submitter_sample_id %in% i]
  sv_tb<-table(df_samp$svtype)
  return(data.table(sample_name=i, project_code, del=sv_tb['DEL'], dup=sv_tb['DUP'], 
        inter=sv_tb['INTER'], inv=sv_tb['INV'], sv_sum=sum(sv_tb)))
}



svtypes_df<-rbindlist(lapply(unique(somatic_subset$sample_name), get_svtypes_pertumor))
svtypes_df[is.na(svtypes_df)]<-0




test_chisq<-function(tt){
  tt_df<-svtypes_df[svtypes_df$project_code==tt]
  tt_df[,not_del:=(sv_sum)-(inv)]
  
  chi_test<-chisq.test(tt_df$inv, tt_df$not_del, correct=F)
  return(data.table(tt, p_val=chi_test$p.value))
}

tt_list<-unique(svtypes_df$project_code)
tt_list<-tt_list[tt_list!=('DLBC-US')]

see_sig_bytt<-rbindlist(lapply(tt_list, test_chisq))


