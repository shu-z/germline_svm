library(data.table)
library(factoextra)




#######cluster misclassified events 

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

