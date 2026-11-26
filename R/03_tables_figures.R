if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','English_United States.utf8'))
suppressPackageStartupMessages({library(gtsummary);library(flextable);library(ggplot2);library(dplyr)})
tabdir<-file.path(output,'tables');dir.create(tabdir,showWarnings=FALSE)
labels<-list(heat5~'Baseline heat per 5 degrees C',age5~'Maternal age per 5 years',smoking~'Smoking yes versus no',deprivation~'Deprivation yes versus no',site~'Illustrative location')
tbl1<-d |>
 select(outcome,age_years,heat14,rain14,smoking,deprivation,site,calendar_year) |>
 tbl_summary(by=outcome,missing='ifany',type=list(calendar_year~'categorical',smoking~'dichotomous',deprivation~'dichotomous'),value=list(smoking~'Yes',deprivation~'Yes'),statistic=list(all_continuous()~'{median} ({p25}, {p75})',all_categorical()~'{n} ({p}%)'),label=list(age_years~'Maternal age years',heat14~'14-day baseline heat degrees C',rain14~'14-day mean rainfall mm',smoking~'Smoking',deprivation~'Deprivation',site~'Illustrative location',calendar_year~'Entry year')) |>
 add_overall(last=FALSE) |> bold_labels()
tbl_uv<-tbl_uvregression(data=d,y=delivery37,method=glm,method.args=list(family=binomial()),include=all_of(covariates),exponentiate=TRUE,tidy_fun=robust_tidy,label=labels,show_single_row=c(smoking,deprivation),hide_n=TRUE,pvalue_fun=label_style_pvalue(digits=3)) |> modify_column_hide(p.value)
tbl_adj<-tbl_regression(primary,exponentiate=TRUE,tidy_fun=robust_tidy,include=all_of(covariates),label=labels,show_single_row=c(smoking,deprivation),pvalue_fun=label_style_pvalue(digits=3))
tbl2<-tbl_merge(list(tbl_uv,tbl_adj),tab_spanner=c('**Unadjusted**','**Adjusted**')) |> bold_labels()
format_p<-function(p)ifelse(p<.001,'<0.001',sprintf('%.3f',p))
sensitivity<-data.frame(Analysis=window_rows$analysis,N=window_rows$n,Events=window_rows$events,`OR (95% CI)`=sprintf('%.2f (%.2f, %.2f)',window_rows$estimate,window_rows$conf.low,window_rows$conf.high),P=format_p(window_rows$p.value),check.names=FALSE)
diagnostic_table<-data.frame(Check=c('Convergence and rank','Separation','Heat linearity on logit scale','Age linearity on logit scale','Collinearity','Grouped goodness of fit','Influential observations'),Finding=c(paste('Converged; rank',primary$rank,'of',length(coef(primary))),if(separation$separation)'Separation detected' else 'No complete or quasi separation detected',paste('Spline versus linear likelihood ratio p',format_p(nonlinearity$p[1])),paste('Spline versus linear likelihood ratio p',format_p(nonlinearity$p[2])),sprintf('Maximum adjusted GVIF %.2f',max(vif$adjusted_GVIF)),sprintf('HL chi-square %.2f on %d df; p %s',hl,hl_df,format_p(hl_p)),paste(sum(cooks>4/nrow(d)),'above Cook distance 4/N; five largest examined')),Interpretation=c('Finite maximum likelihood estimates and estimable design','Linear programming diagnostic; not only a convergence check','Prespecified functional-form sensitivity; no automatic model selection','Prespecified functional-form sensitivity; age need not be normally distributed','Report term degrees of freedom; adjusted GVIF is not the raw one-df VIF','Ten groups; grouping-dependent in-sample diagnostic, not external validation','Primary analysis retains every eligible pregnancy; deletion is sensitivity only'))
flow_table<-flow;names(flow_table)<-c('Preprocessing stage','N')
saveRDS(tbl1,file.path(tabdir,'table1_gtsummary.rds'));saveRDS(tbl2$table_body,file.path(tabdir,'table2_gtsummary_body.rds'))
make_ft<-function(x){
 ft<-if(inherits(x,'gtsummary'))as_flex_table(x) else flextable(x)
 ft<-ft |> font(fontname='Calibri',part='all') |> fontsize(size=9.5,part='all') |> padding(padding=3,part='all') |> bg(bg='#E6EEF2',part='header') |> color(color='black',part='all') |> autofit() |> set_table_properties(layout='autofit',width=1)
 ft<-border_outer(ft,border=officer::fp_border(color='#D9D9D9',width=.5));ft<-border_inner_h(ft,border=officer::fp_border(color='#D9D9D9',width=.5));ft<-border_inner_v(ft,border=officer::fp_border(color='#D9D9D9',width=.5));ft
}
tables<-list(flow=make_ft(flow_table),table1=make_ft(tbl1),table2=make_ft(tbl2),sensitivity=make_ft(sensitivity),diagnostics=make_ft(diagnostic_table))
tables$flow<-paginate(tables$flow,init=TRUE,hdr_ftr=TRUE)
tables$table2<-add_footer_lines(tables$table2,'HC3 robust 95% Wald confidence intervals. Adjusted for age, smoking, deprivation, location, rainfall, annual seasonal terms and calendar year. ORs are associations, not risk ratios or causal effects.')
for(n in names(tables))saveRDS(tables[[n]],file.path(tabdir,paste0(n,'_flextable.rds')))
write_md<-function(x,name){
 x<-if(inherits(x,'gtsummary'))as_tibble(x,col_labels=TRUE) else x
 if(name=='table2')names(x)<-c('Characteristic','Unadjusted OR','95% CI','Adjusted OR','95% CI','Adjusted p')
 cell<-function(v){v<-as.character(v);v[is.na(v)]<-'';v<-gsub('\r?\n','<br>',v);v<-gsub('|','\\|',v,fixed=TRUE);v}
 x[]<-lapply(x,cell);header<-paste0('| ',paste(cell(names(x)),collapse=' | '),' |');rule<-paste0('| ',paste(rep('---',ncol(x)),collapse=' | '),' |')
 body<-apply(x,1,function(r)paste0('| ',paste(r,collapse=' | '),' |'))
 writeLines(c(header,rule,body),file.path(tabdir,paste0(name,'.md')),useBytes=TRUE)
 write.csv(x,file.path(tabdir,paste0(name,'_display.csv')),row.names=FALSE)
}
write_md(tbl1,'table1');write_md(tbl2,'table2');write_md(sensitivity,'sensitivity');write_md(diagnostic_table,'diagnostics');write_md(flow_table,'flow')
packages<-c('gtsummary','flextable','officer','ggplot2','jsonlite','sandwich','car','lpSolve','broom','dplyr','knitr','rmarkdown')
write.csv(data.frame(package=packages,version=sapply(packages,function(x)as.character(packageVersion(x)))),file.path(output,'R_packages.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(output,'R_sessionInfo.txt'))
theme_set(theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face='bold'),plot.subtitle=element_text(size=10)))
save_plot<-function(name,p,width=8,height=4.8){ggsave(file.path(figdir,paste0(name,'.png')),p,width=width,height=height,dpi=300,bg='white');ggsave(file.path(figdir,paste0(name,'.pdf')),p,width=width,height=height,device=cairo_pdf)}
p1<-ggplot(d,aes(site,heat14))+geom_violin(fill='#D6E7EE',colour='#558499',trim=TRUE)+geom_boxplot(width=.15,outlier.shape=NA,fill='white')+coord_flip()+labs(x=NULL,y='Mean daily maximum temperature over 14 days (degrees C)',title='Heat exposure before week-20 entry',subtitle='Public ERA5 weather assigned to 6,000 synthetic pregnancies')
save_plot('01_exposure_distribution',p1,8,3.5)
forest<-primary_terms[primary_terms$term %in% c('heat5','age5','smokingYes','deprivationYes'),];forest$label<-c(heat5='Heat per 5 degrees C',age5='Age per 5 years',smokingYes='Smoking yes vs no',deprivationYes='Deprivation yes vs no')[forest$term]
forest$label<-factor(forest$label,levels=rev(unname(c(heat5='Heat per 5 degrees C',age5='Age per 5 years',smokingYes='Smoking yes vs no',deprivationYes='Deprivation yes vs no'))))
p2<-ggplot(forest,aes(estimate,label))+geom_vline(xintercept=1,linetype=2,colour='grey50')+geom_errorbar(aes(xmin=conf.low,xmax=conf.high),orientation='y',width=.18,colour='#167D9A')+geom_point(colour='#167D9A',size=2.5)+scale_x_log10()+labs(x='Adjusted odds ratio with HC3 95% confidence interval',y=NULL,title='Associations with subsequent delivery before 37 weeks',subtitle='All 6,000 pregnancies retained; covariate-adjusted logistic regression')
save_plot('02_adjusted_associations',p2,8,3.3)
p3<-ggplot(curve,aes(temperature,OR))+geom_ribbon(aes(ymin=lower,ymax=upper),fill='#D6E7EE')+geom_line(colour='#167D9A',linewidth=.9)+geom_hline(yintercept=1,linetype=2)+geom_vline(xintercept=reference,linetype=3)+geom_rug(data=d,aes(x=heat14),inherit.aes=FALSE,alpha=.025)+labs(x='14-day baseline heat (degrees C)',y='Adjusted odds ratio relative to median heat',title='Functional-form assessment of the heat association',subtitle='Natural cubic splines for heat and age; pointwise HC3 95% intervals')
save_plot('03_nonlinear_association',p3)
influence$index<-seq_len(nrow(influence));influence$flag<-influence$id %in% d$id[largest]
p4<-ggplot(influence,aes(index,cooks_distance))+geom_point(colour='#558499',size=.65,alpha=.5)+geom_hline(yintercept=4/nrow(d),linetype=2,colour='#B56528')+geom_point(data=influence[influence$flag,],colour='#B56528',size=1.7)+labs(x='Analysis record',y="Cook's distance",title='Influence diagnostics for the primary association model',subtitle='Orange points: five largest values; all observations remain in the primary fit')
save_plot('04_influence',p4,8,4)
calibration$lower<-with(calibration,((observed_rate+1.96^2/(2*n))-1.96*sqrt(observed_rate*(1-observed_rate)/n+1.96^2/(4*n^2)))/(1+1.96^2/n))
calibration$upper<-with(calibration,((observed_rate+1.96^2/(2*n))+1.96*sqrt(observed_rate*(1-observed_rate)/n+1.96^2/(4*n^2)))/(1+1.96^2/n))
p5<-ggplot(calibration,aes(predicted,observed_rate))+geom_abline(intercept=0,slope=1,linetype=2,colour='grey50')+geom_errorbar(aes(ymin=lower,ymax=upper),width=.002,colour='#167D9A')+geom_point(size=2.5,colour='#167D9A')+coord_equal(xlim=c(.04,.30),ylim=c(.04,.30))+labs(x='Mean fitted event probability within group',y='Observed event proportion',title='Grouped goodness of fit in the analysis cohort',subtitle='Ten groups; descriptive Wilson confidence intervals')
save_plot('05_grouped_fit',p5,6,3.5)
cat('gtsummary tables, native flextables and five ggplot figures saved\n')
