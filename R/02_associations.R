# Prespecified association models and diagnostic checks, not an algorithm contest.
adjustment<-'age5 + smoking + deprivation + site + lograin + season_sin + season_cos + year_c'
full_formula<-as.formula(paste('delivery37 ~ heat5 +',adjustment))
primary<-glm(full_formula,data=d,family=binomial(),control=glm.control(epsilon=1e-10,maxit=100))
crude<-glm(delivery37~heat5,data=d,family=binomial())
maternal<-glm(delivery37~heat5+age5+smoking+deprivation,data=d,family=binomial())
stopifnot(primary$converged,primary$rank==length(coef(primary)),all(is.finite(coef(primary))))
robust_tidy<-function(x,conf.int=TRUE,conf.level=.95,exponentiate=FALSE,...){
 b<-coef(x);v<-sandwich::vcovHC(x,type='HC3');se<-sqrt(diag(v));q<-qnorm(1-(1-conf.level)/2)
 ans<-data.frame(term=names(b),estimate=as.numeric(b),std.error=se,statistic=as.numeric(b/se),p.value=as.numeric(2*pnorm(-abs(b/se))),conf.low=as.numeric(b-q*se),conf.high=as.numeric(b+q*se),row.names=NULL)
 if(exponentiate)ans[,c('estimate','conf.low','conf.high')]<-exp(ans[,c('estimate','conf.low','conf.high')])
 ans
}
primary_terms<-robust_tidy(primary,exponentiate=TRUE);write.csv(primary_terms,file.path(output,'adjusted_odds_ratios.csv'),row.names=FALSE)
covariates<-c('heat5','age5','smoking','deprivation','site')
uv_models<-setNames(lapply(covariates,function(v)glm(reformulate(v,'delivery37'),data=d,family=binomial())),covariates)
uv<-do.call(rbind,lapply(names(uv_models),function(v){a<-robust_tidy(uv_models[[v]],exponentiate=TRUE);a$variable<-v;a}));write.csv(uv,file.path(output,'unadjusted_odds_ratios.csv'),row.names=FALSE)
heat_spline<-update(primary,.~.-heat5+splines::ns(heat5,df=3));age_spline<-update(primary,.~.-age5+splines::ns(age5,df=3));flexible<-update(heat_spline,.~.-age5+splines::ns(age5,df=3))
heat_lrt<-anova(primary,heat_spline,test='Chisq');age_lrt<-anova(primary,age_spline,test='Chisq')
nonlinearity<-data.frame(variable=c('Heat','Maternal age'),df=c(heat_lrt$Df[2],age_lrt$Df[2]),LR=c(heat_lrt$Deviance[2],age_lrt$Deviance[2]),p=c(heat_lrt[['Pr(>Chi)']][2],age_lrt[['Pr(>Chi)']][2]))
write.csv(nonlinearity,file.path(output,'nonlinearity_tests.csv'),row.names=FALSE)
windows<-list('14-day exposure; unadjusted'=crude,'14-day exposure; maternal adjusted'=maternal,'14-day exposure; fully adjusted'=primary,'7-day exposure; fully adjusted'=update(primary,.~.-heat5+heat7_5),'Days 8-14; fully adjusted'=update(primary,.~.-heat5+previous7_5),'Site-specific seasonality'=update(primary,.~.+site:season_sin+site:season_cos),'Exploratory delivery before 34 weeks'=update(primary,delivery34~.))
cooks<-cooks.distance(primary);hat<-hatvalues(primary);largest<-order(cooks,decreasing=TRUE)[1:5]
windows[['Omit five highest Cook distances']]<-glm(full_formula,data=d[-largest,],family=binomial())
window_rows<-do.call(rbind,lapply(names(windows),function(name){m<-windows[[name]];term<-if(grepl('^7-day',name))'heat7_5' else if(grepl('^Days 8',name))'previous7_5' else 'heat5';a<-robust_tidy(m,exponentiate=TRUE);a<-a[a$term==term,,drop=FALSE];a$analysis<-name;a$n<-nobs(m);a$events<-sum(model.response(model.frame(m)));a}))
write.csv(window_rows,file.path(output,'exposure_associations.csv'),row.names=FALSE)
vif<-car::vif(primary);vif<-data.frame(term=rownames(vif),GVIF=vif[,1],df=vif[,2],adjusted_GVIF=vif[,3],row.names=NULL);write.csv(vif,file.path(output,'collinearity_gvif.csv'),row.names=FALSE)
# Linear programming: a nonzero direction with nonnegative transformed margins
# and at least one strict margin characterises quasi or complete separation.
separation_check<-function(x,y){
 A<-x*(2*y-1);B<-cbind(A,-A);p<-ncol(B)
 quasi<-lpSolve::lp('min',rep(0,p),rbind(B,colSums(B)),c(rep('>=',nrow(B)),'>='),c(rep(0,nrow(B)),1))
 complete<-lpSolve::lp('min',rep(0,p),B,rep('>=',nrow(B)),rep(1,nrow(B)))
 list(separation=quasi$status==0,complete=complete$status==0,quasi_status=quasi$status,complete_status=complete$status)
}
separation<-separation_check(model.matrix(primary),d$delivery37)
pred<-fitted(primary);bins<-cut(pred,breaks=unique(quantile(pred,seq(0,1,.1))),include.lowest=TRUE)
calibration<-do.call(rbind,lapply(split(seq_along(pred),bins),function(j)data.frame(n=length(j),observed=sum(d$delivery37[j]),expected=sum(pred[j]),predicted=mean(pred[j]),observed_rate=mean(d$delivery37[j]),mean_residual=mean(d$delivery37[j]-pred[j]))))
hl<-sum((calibration$observed-calibration$expected)^2/(calibration$expected*(1-calibration$expected/calibration$n)))
hl_df<-nrow(calibration)-2;hl_p<-pchisq(hl,hl_df,lower.tail=FALSE)
write.csv(calibration,file.path(output,'grouped_fit_check.csv'),row.names=FALSE)
influence<-data.frame(id=d$id,cooks_distance=cooks,leverage=hat,deviance_residual=residuals(primary,type='deviance'),pearson_residual=residuals(primary,type='pearson'),fitted_probability=pred)
write.csv(influence,file.path(output,'influence_diagnostics.csv'),row.names=FALSE)
grid<-seq(quantile(d$heat14,.01),quantile(d$heat14,.99),length.out=120);reference<-median(d$heat14)
new<-d[rep(1,length(grid)),];new$heat5<-grid/5;ref<-new;ref$heat5<-reference/5
mx<-model.matrix(delete.response(terms(flexible)),new,contrasts.arg=flexible$contrasts);mr<-model.matrix(delete.response(terms(flexible)),ref,contrasts.arg=flexible$contrasts);delta<-mx-mr;v<-sandwich::vcovHC(flexible,type='HC3');b<-coef(flexible)
eta<-drop(delta%*%b);se<-sqrt(rowSums((delta%*%v)*delta))
curve<-data.frame(temperature=grid,OR=exp(eta),lower=exp(eta-1.96*se),upper=exp(eta+1.96*se),reference=reference)
write.csv(curve,file.path(output,'nonlinear_heat_association.csv'),row.names=FALSE)
iqr_data<-d[c(1,1),];iqr_data$heat5<-as.numeric(quantile(d$heat14,c(.25,.75)))/5
iqr_matrix<-model.matrix(delete.response(terms(flexible)),iqr_data,contrasts.arg=flexible$contrasts);con<-iqr_matrix[2,]-iqr_matrix[1,];iqr_eta<-sum(con*b);iqr_se<-sqrt(drop(t(con)%*%v%*%con))
diagnostics<-list(n=nrow(d),events=sum(d$delivery37),converged=primary$converged,model_rank=primary$rank,coefficients=length(coef(primary)),full_rank=primary$rank==length(coef(primary)),min_fitted=min(pred),max_fitted=max(pred),separation=separation,heat_nonlinearity_p=nonlinearity$p[1],age_nonlinearity_p=nonlinearity$p[2],max_adjusted_gvif=max(vif$adjusted_GVIF),HL_statistic=hl,HL_df=hl_df,HL_p=hl_p,cook_threshold=4/nrow(d),above_cook_threshold=sum(cooks>4/nrow(d)),max_cook=max(cooks),leverage_threshold=2*primary$rank/nrow(d),above_leverage_threshold=sum(hat>2*primary$rank/nrow(d)),influence_deleted_ids=d$id[largest],residual_deviance=deviance(primary),residual_df=df.residual(primary),inference='HC3 robust Wald intervals; no posterior quantities',spline_iqr=list(low=as.numeric(quantile(d$heat14,.25)),high=as.numeric(quantile(d$heat14,.75)),OR=exp(iqr_eta),lower=exp(iqr_eta-1.96*iqr_se),upper=exp(iqr_eta+1.96*iqr_se)))
jsonlite::write_json(diagnostics,file.path(output,'model_diagnostics.json'),auto_unbox=TRUE,pretty=TRUE,digits=12)
saveRDS(list(primary=primary,crude=crude,maternal=maternal,uv=uv_models,heat_spline=heat_spline,age_spline=age_spline,flexible=flexible,windows=windows),file.path(output,'association_models.rds'))
cat('Association models complete; adjusted heat OR',primary_terms$estimate[primary_terms$term=='heat5'],'\n')
