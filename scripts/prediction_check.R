args<-commandArgs(trailingOnly=TRUE);root<-args[1]
d<-read.csv(file.path(root,'data/landmark_design.csv'));x<-as.matrix(d[,grepl('^x[0-9]+$',names(d))]);s<-d$train==1
f<-glm.fit(x[s,],d$outcome[s],family=binomial(),control=glm.control(epsilon=1e-10,maxit=100))
write.csv(data.frame(term=0:(ncol(x)-1),coefficient=f$coefficients),file.path(root,'results/r_landmark_logistic.csv'),row.names=FALSE)
