args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[1] else '.',winslash='/')
manifest <- if(length(args)>1) args[2] else 'data/simulation_manifest.csv'
output <- if(length(args)>2) args[3] else 'results/r_fits.csv'
w <- read.csv(file.path(root,'data/weather.csv'))
weather <- lapply(0:5,function(s) {d<-w[w$site==s,];list(date=as.Date(d$date),ct=c(0,cumsum(d$tmax_c)),cr=c(0,cumsum(d$rain_mm)))})
heat_basis <- function(t0,t1) {a<-pmax(t0-25,0)/5;b<-pmax(t1-25,0)/5;cbind(a,a*a,b,b*b)}
history <- function(c) {
 p<-rep(seq_len(nrow(c)),each=17);k<-rep(0:16,nrow(c));day<-c$entry_day[p]+7*k
 t0<-t1<-rain<-doy<-year<-numeric(length(p))
 for(s in 0:5){sel<-c$site[p]==s;d<-day[sel];z<-weather[[s+1]]
  t0[sel]<-(z$ct[d+1]-z$ct[d-7+1])/7;t1[sel]<-(z$ct[d-7+1]-z$ct[d-14+1])/7
  rain[sel]<-(z$cr[d+1]-z$cr[d-7+1])/7;date<-z$date[d+1];doy[sel]<-as.integer(format(date,'%j'));year[sel]<-as.integer(format(date,'%Y'))
 }
 ga<-(k+20-28)/4
 base<-cbind(1,ga,ga*ga,c$age_c[p],c$smoking[p],c$deprivation[p],log1p(rain),sin(2*pi*doy/365.25),cos(2*pi*doy/365.25),year-2022,sapply(1:5,function(s)as.numeric(c$site[p]==s)))
 list(p=p,k=k,base=base,t0=t0,t1=t1)
}
mf<-read.csv(file.path(root,manifest));records<-list();coeff<-list();counter<-0;start<-proc.time()[3]
for(i in seq_len(nrow(mf))){r<-mf[i,];c<-read.csv(gzfile(file.path(root,r$path)));h<-history(c);p<-h$p;week<-h$k+20;y<-as.numeric(week==c$birth_week[p]);cached<-NULL
 for(model in c('nonlinear','linear','oracle')){
  oracle<-model=='oracle';keep<-week<=pmin(c$birth_week[p],36)
  if(!oracle)keep<-keep & bitwAnd(c$missing_mask[p],bitwShiftL(1L,h$k))==0
  t0<-h$t0+if(oracle)c$personal_offset[p] else 0;t1<-h$t1+if(oracle)c$personal_offset[p] else 0
  heat<-if(model=='linear')matrix((t0-25)/5,ncol=1) else heat_basis(t0,t1);x<-cbind(h$base,heat)[keep,,drop=FALSE];yy<-y[keep];counter<-counter+1
  result<-tryCatch({
   f<-if(model=='oracle' && r$scenario %in% c('null','moderate','nonlinear'))cached else glm.fit(x,yy,family=binomial('cloglog'),control=glm.control(epsilon=1e-10,maxit=100));if(model=='nonlinear')cached<-f
   v<-summary.glm(f)$cov.unscaled
   contrast<-c(rep(0,15),rep(1,ncol(x)-15));theta<-sum(contrast*f$coefficients);se<-sqrt(drop(t(contrast)%*%v%*%contrast))
   if(r$replicate==1){for(j in seq_len(ncol(x)))coeff[[length(coeff)+1]]<-data.frame(scenario=r$scenario,model=model,term=j-1,estimate=f$coefficients[j],se=sqrt(diag(v))[j])}
   data.frame(replicate=r$replicate,scenario=r$scenario,model=model,n=r$n,events=r$events,rows=length(yy),observed_events=sum(yy),truth=r$truth,estimate=theta,se=se,lower=theta-1.96*se,upper=theta+1.96*se,p=2*pnorm(-abs(theta/se)),converged=f$converged,loglik=sum(dbinom(yy,1,f$fitted.values,log=TRUE)),error='')
  },error=function(e)data.frame(replicate=r$replicate,scenario=r$scenario,model=model,n=r$n,events=r$events,rows=length(yy),observed_events=sum(yy),truth=r$truth,estimate=NA,se=NA,lower=NA,upper=NA,p=NA,converged=FALSE,loglik=NA,error=conditionMessage(e)))
  records[[counter]]<-result
 }
 if(r$replicate%%10==0 && r$scenario==tail(mf$scenario,1))cat('R replicate',r$replicate,'seconds',round(proc.time()[3]-start),'\n')
}
dir.create(file.path(root,'results'),showWarnings=FALSE)
write.csv(do.call(rbind,records),file.path(root,output),row.names=FALSE)
write.csv(do.call(rbind,coeff),sub('.csv','_coefficients.csv',file.path(root,output),fixed=TRUE),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(root,'provenance/R_sessionInfo.txt'))
cat('R fits',counter,'seconds',round(proc.time()[3]-start),'\n')
