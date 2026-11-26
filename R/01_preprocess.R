# Read, audit and derive the association cohort from the saved source files.
config<-jsonlite::fromJSON(file.path(root,'association_config.json'))
output<-file.path(root,'results/association');figdir<-file.path(root,'figures/association')
dir.create(output,recursive=TRUE,showWarnings=FALSE);dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
raw<-read.csv(file.path(root,config$source_cohort),check.names=FALSE)
weather<-read.csv(file.path(root,config$source_weather),check.names=FALSE)
stopifnot(!anyDuplicated(raw$id),!anyDuplicated(weather[c('site','date')]),nrow(raw)==6000,nrow(weather)==13152)
weather$date<-as.Date(weather$date)
stopifnot(all(is.finite(weather$tmax_c)),all(is.finite(weather$tmin_c)),all(is.finite(weather$rain_mm)),all(weather$tmin_c<=weather$tmax_c),all(weather$rain_mm>=0),all(weather$tmin_c>=-30),all(weather$tmax_c<=55))
for(s in 0:5){z<-weather[weather$site==s,];stopifnot(nrow(z)==2192,all(diff(z$date)==1))}
stopifnot(all(raw$site %in% 0:5),all(raw$birth_week %in% 20:37),all(raw$smoking %in% 0:1),all(raw$deprivation %in% 0:1),all(raw$missing_mask==0),all(raw$personal_offset==0))
d<-data.frame(id=raw$id,site_code=raw$site,entry_date=as.Date('2019-01-01')+raw$entry_day,age_years=30+10*raw$age_c,smoking=factor(raw$smoking,levels=0:1,labels=c('No','Yes')),deprivation=factor(raw$deprivation,levels=0:1,labels=c('No','Yes')),birth_week=raw$birth_week)
stopifnot(all(d$age_years>=18 & d$age_years<=44))
d$heat7<-d$previous7<-d$rain14<-NA_real_
for(s in 0:5){
 z<-weather[weather$site==s,];ct<-c(0,cumsum(z$tmax_c));cr<-c(0,cumsum(z$rain_mm));sel<-which(d$site_code==s);k<-match(d$entry_date[sel],z$date)
 stopifnot(!anyNA(k),all(k>14))
 # ct[k] sums the days strictly before entry; the entry day is excluded.
 d$heat7[sel]<-(ct[k]-ct[k-7])/7;d$previous7[sel]<-(ct[k-7]-ct[k-14])/7
 d$rain14[sel]<-(cr[k]-cr[k-14])/14
}
d$heat14<-(d$heat7+d$previous7)/2;d$heat5<-d$heat14/5;d$heat7_5<-d$heat7/5;d$previous7_5<-d$previous7/5;d$age5<-(d$age_years-30)/5
d$lograin<-log1p(d$rain14);d$calendar_year<-as.integer(format(d$entry_date,'%Y'));d$year_c<-d$calendar_year-2022
d$doy<-as.integer(format(d$entry_date,'%j'));d$season_sin<-sin(2*pi*d$doy/365.25);d$season_cos<-cos(2*pi*d$doy/365.25)
site_names<-c('Johannesburg','Durban','Cape Town','Harare','Nairobi','Abidjan')
d$site<-factor(site_names[d$site_code+1],levels=site_names)
d$delivery37<-as.integer(d$birth_week<37);d$delivery34<-as.integer(d$birth_week<34)
d$outcome<-factor(d$delivery37,levels=0:1,labels=c('No delivery before 37 weeks','Delivery before 37 weeks'))
d$exposure_start<-d$entry_date-14;d$exposure_end<-d$entry_date-1
d$event_interval_start<-as.Date(ifelse(d$delivery37==1,as.numeric(d$entry_date+7*(d$birth_week-20)),NA),origin='1970-01-01')
stopifnot(all(d$exposure_end<d$entry_date),all(d$event_interval_start[d$delivery37==1]>=d$entry_date[d$delivery37==1]))
variables<-c('delivery37','heat5','age5','smoking','deprivation','site','lograin','season_sin','season_cos','year_c')
missingness<-data.frame(variable=names(d),missing=sapply(d,function(x)sum(is.na(x))),total=nrow(d))
missingness$percent<-100*missingness$missing/missingness$total
write.csv(missingness,file.path(output,'missingness_audit.csv'),row.names=FALSE)
keep<-complete.cases(d[,variables]);flow<-data.frame(stage=c('Source synthetic pregnancies','Duplicate identifiers','Outside prespecified age or event ranges','Incomplete primary model records','Included in primary analysis','Delivery before week 37','Delivery before week 34'),n=c(nrow(raw),sum(duplicated(raw$id)),0,sum(!keep),sum(keep),sum(d$delivery37[keep]),sum(d$delivery34[keep])))
d<-d[keep,];stopifnot(nrow(d)==6000,sum(d$delivery37)==796)
write.csv(flow,file.path(output,'cohort_flow.csv'),row.names=FALSE);write.csv(d,file.path(output,'analysis_cohort.csv'),row.names=FALSE)
saveRDS(d,file.path(output,'analysis_cohort.rds'))
dictionary<-data.frame(variable=c('delivery37','delivery34','birth_week','heat14','heat7','previous7','heat5','age_years','age5','smoking','deprivation','site','rain14','season_sin','season_cos','year_c'),definition=c('Primary generated delivery event: 0 none through week 36; 1 event in weeks 20-36','Exploratory generated delivery event before week 34','Generated interval; 37 denotes no earlier event and is not an exact term delivery week','Mean daily maximum temperature during entry days -14 through -1 in degrees C','Temperature on entry days -7 through -1','Temperature on entry days -14 through -8','Primary exposure divided by 5; coefficient represents a 5 degree C contrast','Reconstructed as 30 + 10 times age_c','Age centred at 30 divided by 5','Synthetic binary smoking status','Synthetic binary deprivation status','Fixed illustrative location; Johannesburg reference','Mean precipitation during entry days -14 through -1 in mm','Annual sine at week-20 entry','Annual cosine at week-20 entry','Entry year centred at 2022'))
write.csv(dictionary,file.path(output,'data_dictionary.csv'),row.names=FALSE)
preprocess<-list(source_n=nrow(raw),analysis_n=nrow(d),events37=sum(d$delivery37),events34=sum(d$delivery34),primary_missing_records=sum(!keep),duplicate_ids=sum(duplicated(raw$id)),weather_rows=nrow(weather),exposure_days=14,heat_median=median(d$heat14),heat_q25=unname(quantile(d$heat14,.25)),heat_q75=unname(quantile(d$heat14,.75)),entry_min=as.character(min(d$entry_date)),entry_max=as.character(max(d$entry_date)),site_labels_are='Illustrative weather locations, not clinical recruitment sites',outcome_limitation='Generator does not distinguish live birth from pregnancy loss; do not interpret as a validated obstetric preterm-birth endpoint')
jsonlite::write_json(preprocess,file.path(output,'preprocessing_summary.json'),auto_unbox=TRUE,pretty=TRUE)
cat('Preprocessing complete:',nrow(d),'pregnancies and',sum(d$delivery37),'primary events\n')
