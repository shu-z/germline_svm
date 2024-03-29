library(data.table)
library(e1071)
library(caTools)
library(ROCR)
library(stringr)
library(ggfortify)
library(corrplot)
library(Hmisc)


#use after train/test sets made from svm_main 



#separate continuous and categorical 
continuous_vars<-c('homlen', 'insertion_len', 'SPAN', 'gnomad_dist', 'hom_gc', 
                   'insertion_gc', 'line_dist', 'sine_dist', 'num_sv_sample')

binary_vars<-c('del', 'dup', 'inv', 'inter', 'CN_annot', 'exon_annot')


#####################################
#look at correlations of continuous features 
somatic_continuous<-somatic_subset[, ..continuous_vars]
germline_continuous<-germline_subset[, ..continuous_vars]

#use rcorr to get the pvalues
cor_somatic <-rcorr(as.matrix(somatic_continuous), type='spearman')
somatic_r<-(cor_somatic$r)
#View(cor_somatic$P)
cor_germline <- rcorr(as.matrix(germline_continuous), type='spearman')
germline_r<-(cor_germline$r)
#View(cor_germline$P)
#cor_all <- rcorr(as.matrix(df[, ..features_toscale]), type='spearman')
#View(cor_all$P)


#pdf('/Users/shu/germline_svm/figs/20220426_newannot_corrplot_pearson_continuousfeatures.pdf', width=8, height=8)
corrplot(cor_somatic, type='upper', addCoef.col = 'grey50', number.cex = 0.75, mar=c(0,0,2,0)) 
mtext("somatic, pearson, continuous features", line=2, cex=1)
corrplot(cor_germline, type='upper', addCoef.col = 'grey50', number.cex = 0.75, mar=c(0,0,2,0))
mtext("germline, pearson, continuous features", line=2, cex=1)
#dev.off()

#look at all
#pdf('/Users/shu/germline_svm/figs/20220426_newannot_corrplotcombined_pearson_allfeatures.pdf', width=8, height=8)
corrplot(cor_all, type='upper', addCoef.col = 'grey50', number.cex = 0.75, mar=c(0,0,2,0)) 
mtext("10k germline, 10k somatic, pearson, all features", line=2, cex=1)
#dev.off()



#transform correlation rho into Z scores
r_to_z<-function(line){
  
  #idx is indices of feature from continuous_vars
  idx<-upper_indices[line,]
  
  r_g<-(germline_r[idx$row,idx$col])
  germline_z<-0.5*(log(1+r_g) - log(1-r_g))
  
  r_so<-(somatic_r[idx$row,idx$col])
  somatic_z<-0.5*(log(1+r_so) - log(1-r_so))
  
  #combine germline and somatic Z
  z_divide<-1/(nrow(somatic_subset)-3) + 1/(nrow(germline_subset)-3)
  z_test<-(somatic_z-germline_z)/(sqrt(z_divide))
  
  #get one-sided pval based on sign of Z
  ifelse(z_test>0, pval<-pnorm(z_test, lower.tail=F), pval<-pnorm(z_test, lower.tail=T))
  
  return(data.table(continuous_vars[idx$row], continuous_vars[idx$col], 
                    germline_z, somatic_z, z_divide, pval))
}


#get all upper triangular indices (all possible feature pairs), and run rho-to-Z
g <- expand.grid(row = 1:nrow(germline_r), col = 1:ncol(germline_r)) # grid
upper_indices<-g[upper.tri(germline_r, diag = F), ]
z_scores<-rbindlist(lapply(1:nrow(upper_indices), r_to_z))

###############################################


######## compare binary vs binary 

somatic_binary<-somatic_subset[, ..binary_vars]
germline_binary<-germline_subset[, ..binary_vars]

#chi sq test for binary vs binary 
binary_chisq<-function(line, df){
  idx<-upper_indices_binary[line,]
  i<-as.numeric(idx$row)
  j<-as.numeric(idx$col)

  test<-chisq.test(unlist(df[, ..i]), unlist(df[, ..j]))
  return(data.table(binary_vars[i], binary_vars[j], pval=test$p.value))
}

#get all upper triangular indices (all possible feature pairs), and run chi sq
g <- expand.grid(row = 1:ncol(somatic_binary), col = 1:ncol(somatic_binary)) 
empty_mat<-matrix(NA, ncol(somatic_binary), ncol(somatic_binary))
upper_indices_binary<-g[upper.tri(empty_mat, diag = F), ]

germline_binary_pval<-rbindlist(lapply(1:nrow(upper_indices_binary), binary_chisq, df=germline_binary))
somatic_binary_pval<-rbindlist(lapply(1:nrow(upper_indices_binary), binary_chisq, df=somatic_binary))
germline_binary_pval[,qval:=p.adjust(pval, method='fdr')]
somatic_binary_pval[,qval:=p.adjust(pval, method='fdr')]

###############################################


######### mann whitney for binary vs. continous vars 

binary_mw<-function(line, df_binary, df_continuous){
  idx<-binary_cont_indices[line,]
  #i is the binary feature, j is continuous 
  i<-as.numeric(idx$row)
  j<-as.numeric(idx$col)
  
  df_pair<-cbind(df_binary[, ..i], df_continuous[, ..j], fill=T)
  
  #get indices where binary feature = 0
  binary_0<-as.vector(df_pair[,1]==0)
  #two dfs for binary=1 and binary=0
  df_pair_0<-df_pair[binary_0,]
  df_pair_1<-df_pair[(!binary_0),]
  
  #mann whitney test
  test<-wilcox.test(unlist(df_pair_0[,2]), unlist(df_pair_1[,2]))
  return(data.table(binary_vars[i], continuous_vars[j], pval=test$p.value))
}

#get all upper triangular indices, and run 
binary_cont_indices <- expand.grid(row = 1:ncol(somatic_binary), col = 1:ncol(somatic_continuous)) 

germline_mw_pval<-rbindlist(lapply(1:nrow(binary_cont_indices), binary_mw, df_binary=germline_binary, df_continuous=germline_continuous))
somatic_mw_pval<-rbindlist(lapply(1:nrow(binary_cont_indices), binary_mw, df=somatic_binary, df_continuous=somatic_continuous))
germline_mw_pval[,qval:=p.adjust(pval, method='fdr')]
somatic_mw_pval[,qval:=p.adjust(pval, method='fdr')]





