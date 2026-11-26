args<-commandArgs(trailingOnly=TRUE);root<-args[1];packages<-c('ggplot2','ranger','gbm','glmnet','pROC','jsonlite','rmarkdown','knitr')
write.csv(data.frame(package=packages,version=sapply(packages,function(p)as.character(packageVersion(p)))),file.path(root,'provenance/R_packages.csv'),row.names=FALSE)
