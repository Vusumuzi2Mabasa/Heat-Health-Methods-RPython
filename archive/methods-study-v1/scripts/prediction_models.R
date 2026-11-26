args<-commandArgs(trailingOnly=TRUE);root<-args[1]
suppressPackageStartupMessages({library(ranger);library(gbm);library(glmnet);library(pROC)})
d<-read.csv(file.path(root,'data/landmark_design.csv'));c<-read.csv(file.path(root,'data/illustrative_cohort.csv'))
x<-as.matrix(d[,paste0('x',1:9)]);y<-d$outcome;date<-as.Date('2019-01-01')+c$entry_day;year<-as.integer(format(date,'%Y'));eligible<-c$birth_week>=28
early<-date+119<=as.Date('2022-12-31');late<-date>=as.Date('2023-01-01')
split<-list(train=eligible&c$site<4&early,temporal=eligible&c$site<4&late,geographic=eligible&c$site>=4&early,combined_shift=eligible&c$site>=4&late);tr<-split$train
set.seed(20261126);train<-data.frame(y=factor(y[tr],levels=0:1),x[tr,]);colnames(train)<-c('y',paste0('f',1:9));all<-data.frame(x);names(all)<-paste0('f',1:9)
ridge<-glmnet(x[tr,],y[tr],family='binomial',alpha=0,lambda=1/sum(tr),standardize=TRUE,thresh=1e-10)
maternal<-glmnet(x[tr,1:3],y[tr],family='binomial',alpha=0,lambda=1/sum(tr),standardize=TRUE,thresh=1e-10)
rf<-ranger(y~.,data=train,num.trees=400,min.node.size=25,min.bucket=25,mtry=6,probability=TRUE,num.threads=1,seed=20261126)
boost<-gbm.fit(x[tr,],y[tr],distribution='bernoulli',n.trees=150,interaction.depth=2,shrinkage=.03,n.minobsinnode=25,bag.fraction=1,verbose=FALSE)
pred<-list('Maternal baseline'=as.numeric(predict(maternal,x[,1:3],type='response',s=1/sum(tr))),Logistic=as.numeric(predict(ridge,x,type='response',s=1/sum(tr))),'Random forest'=predict(rf,all)$predictions[,'1'],'Gradient boosting'=as.numeric(predict(boost,x,n.trees=150,type='response')))
metrics<-list();rocrows<-list();calrows<-list();predrows<-list();count<-0
for(model in names(pred))for(test in names(split)[-1]){
 sel<-split[[test]];yy<-y[sel];pr<-pred[[model]][sel];z<-qlogis(pmin(pmax(pr,1e-6),1-1e-6));f<-glm(yy~z,family=binomial());inter<-glm(yy~1+offset(z),family=binomial());roc<-pROC::roc(yy,pr,direction='<',quiet=TRUE)
 # Average precision integrates precision over recall increments at distinct tied scores.
 scores<-sort(unique(pr),decreasing=TRUE);tp<-sapply(scores,function(t)sum(yy[pr>=t]));fp<-sapply(scores,function(t)sum(yy[pr>=t]==0));rec<-tp/sum(yy);ap<-sum(diff(c(0,rec))*tp/(tp+fp))
 aucs<-replicate(300,{j<-sample(seq_along(yy),replace=TRUE);if(length(unique(yy[j]))==2)as.numeric(pROC::auc(yy[j],pr[j],direction='<',quiet=TRUE)) else NA_real_});ci<-quantile(aucs,c(.025,.975),na.rm=TRUE)
 count<-count+1;metrics[[count]]<-data.frame(model=model,test=test,n=length(yy),events=sum(yy),prevalence=mean(yy),auc=as.numeric(pROC::auc(roc)),auc_lower=ci[1],auc_upper=ci[2],average_precision=ap,brier=mean((yy-pr)^2),calibration_intercept=coef(inter)[1],calibration_slope=coef(f)[2],mean_prediction=mean(pr))
 predrows[[count]]<-data.frame(id=c$id[sel],model=model,test=test,outcome=yy,prediction=pr)
 rocrows[[count]]<-data.frame(model=model,test=test,fpr=1-roc$specificities,tpr=roc$sensitivities)
 bins<-cut(rank(pr,ties.method='first'),breaks=seq(0,length(pr),length.out=6),include.lowest=TRUE)
 calrows[[count]]<-do.call(rbind,lapply(split(seq_along(yy),bins),function(j)data.frame(model=model,test=test,mean_prediction=mean(pr[j]),observed=mean(yy[j]),n=length(j))))
}
write.csv(do.call(rbind,metrics),file.path(root,'results/r_prediction_metrics.csv'),row.names=FALSE)
write.csv(do.call(rbind,rocrows),file.path(root,'results/r_roc.csv'),row.names=FALSE)
write.csv(do.call(rbind,calrows),file.path(root,'results/r_calibration.csv'),row.names=FALSE)
write.csv(do.call(rbind,predrows),file.path(root,'results/r_landmark_predictions.csv'),row.names=FALSE)
saveRDS(list(maternal=maternal,ridge=ridge,random_forest=rf,gradient_boosting=boost,features=colnames(all),train_ids=c$id[tr],seed=20261126),file.path(root,'results/frozen_R_landmark_models.rds'))
print(do.call(rbind,metrics))
