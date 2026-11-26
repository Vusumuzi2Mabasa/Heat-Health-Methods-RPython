# Both manuscript outputs use the same R-derived narrative and table objects.
output<-file.path(root,'results/association');tabdir<-file.path(output,'tables')
d<-readRDS(file.path(output,'analysis_cohort.rds'))
primary_terms<-read.csv(file.path(output,'adjusted_odds_ratios.csv'))
window_rows<-read.csv(file.path(output,'exposure_associations.csv'))
nonlinearity<-read.csv(file.path(output,'nonlinearity_tests.csv'));vif<-read.csv(file.path(output,'collinearity_gvif.csv'))
cooks<-read.csv(file.path(output,'influence_diagnostics.csv'))$cooks_distance
format_p<-function(p)ifelse(p<.001,'<0.001',sprintf('%.3f',p))
p_text<-function(p)ifelse(p<.001,'p < 0.001',sprintf('p = %.3f',p))
prep<-jsonlite::fromJSON(file.path(output,'preprocessing_summary.json'));diag<-jsonlite::fromJSON(file.path(output,'model_diagnostics.json'))
report_rank<-diag$model_rank;hl<-diag$HL_statistic;hl_df<-diag$HL_df;hl_p<-diag$HL_p
effect_text<-function(term){z<-primary_terms[primary_terms$term==term,];sprintf('%.2f (95%% CI %.2f-%.2f; %s)',z$estimate,z$conf.low,z$conf.high,p_text(z$p.value))}
heat<-effect_text('heat5');age<-effect_text('age5');smoke<-effect_text('smokingYes');depr<-effect_text('deprivationYes')
raw_heat<-window_rows[window_rows$analysis=='14-day exposure; unadjusted',]
influence_row<-window_rows[window_rows$analysis=='Omit five highest Cook distances',]
season_row<-window_rows[window_rows$analysis=='Site-specific seasonality',]
early_row<-window_rows[window_rows$analysis=='Exploratory delivery before 34 weeks',]
body<-r"(## Abstract

**Background.** Maternal heat exposure may be associated with pregnancy outcomes, but an association analysis requires a fixed exposure window, clearly defined endpoints and appropriate covariate adjustment. This study investigated early gestational heat and subsequent delivery in a synthetic cohort linked to public historical weather.

**Methods.** We analysed 6,000 synthetic pregnancies ongoing at gestational week 20. The primary exposure was mean daily maximum temperature during the 14 completed days before entry. The primary outcome was generated delivery during weeks 20-36, coded 1, versus no delivery event through week 36, coded 0. Data preprocessing, descriptive analysis, logistic regression and diagnostics were conducted in R. Adjusted odds ratios used HC3 sandwich standard errors. Tables were generated with gtsummary and rendered in Word using flextable.

@@ABSTRACT_RESULTS@@

**Conclusion.** The crude positive heat association was not retained after full covariate adjustment in the primary linear model. Age, smoking and deprivation showed positive adjusted associations in this generated cohort. These results demonstrate an applied association analysis of simulated health data; they do not establish clinical effects or a protective effect of heat.

## Introduction

Heat exposure is relevant to maternal and newborn health research. A systematic review reported associations between higher temperatures during pregnancy and outcomes including preterm birth [1]. Such evidence motivates careful exposure and endpoint definitions. A temperature measurement obtained before an observation period answers a different question from a whole-pregnancy average or a time-varying exposure series.

The present study uses a single synthetic cohort to examine an epidemiological association from data preparation through interpretation. Public ERA5 weather provides the environmental measurements, while all maternal characteristics and delivery outcomes are simulated. The analysis focuses on clinical questions, descriptive findings and covariate-adjusted estimates. The original statistical-methods experiment is preserved as a separate archived release and is not the subject of this manuscript.

### Aim and objectives

The aim was to investigate the association between heat exposure before week-20 entry and subsequent delivery before week 37 in a simulated pregnancy cohort.

The objectives were to describe maternal characteristics, environmental exposure and delivery events; estimate unadjusted and adjusted heat associations; assess functional-form and other model assumptions; and examine prespecified exposure-window and adjustment sensitivities. Maternal age, smoking and deprivation associations were reported as secondary covariate findings. Delivery before week 34 was an exploratory secondary endpoint.

## Methods

### Study design and data sources

This was an association analysis of 6,000 synthetic pregnancies assigned equally across six illustrative weather locations: Johannesburg, Durban, Cape Town, Harare, Nairobi and Abidjan. The locations are environmental reference points, not clinical recruitment sites. The cohort is not a representative sample of any country or maternity service.

The existing source cohort was retained without regenerating outcomes. Its delivery intervals had been generated using gestational time, maternal characteristics and evolving weather exposure. Consequently, the baseline association evaluated here is not the same estimand as the generator's weekly heat coefficient. No parameters or participant records were changed to make the associations statistically significant.

Daily maximum temperature, minimum temperature and precipitation were obtained previously from ERA5 through the Open-Meteo Historical Weather API [2]. The saved panel covered 2019-2024 with 2,192 consecutive days per location and 13,152 site-days. Reanalysis represents gridded environmental estimates, not personal heat exposure. Health data and identifiers are entirely synthetic.

### Preprocessing and quality assessment

R imported the source CSV files and verified unique participant identifiers and site-date keys. Checks covered valid binary covariates, age and generated delivery-week ranges, complete daily weather histories, finite measurements, ordered minimum and maximum temperatures, plausible temperature bounds and nonnegative precipitation. Source age was reconstructed in years as 30 plus ten times the stored centred-age variable. Binary smoking and deprivation were encoded as labelled factors, with No as reference; outcome coding remained 0 and 1.

Exposure windows were linked by location and date. Direct summation of daily values in 60 distributed records independently checked the cumulative-sum implementation. Analysis-variable missingness was audited before selecting complete records. No duplicate, invalid or incomplete primary records were identified; no imputation or automatic outlier removal was performed. Event-interval start dates are structurally missing for non-events and were not predictors or complete-case selection variables.

@@FLOW_TABLE@@

The final two audit rows are outcome counts, not additional exclusions. The analysis dataset, data dictionary, missingness audit and cohort accounting are retained with the R outputs. All original result files and both source datasets were checked against the frozen archive.

### Exposure and outcome definitions

Let $d_i$ be the date of week-20 entry and $T_{s_i}(d)$ the daily maximum temperature at the assigned location. The primary exposure was

$$
H_i=\frac{1}{14}\sum_{j=1}^{14}T_{s_i}(d_i-j).
$$

The entry day was excluded, and every exposure date preceded every possible generated delivery interval. This fixed window avoids defining exposure duration using the subsequent delivery date. Two seven-day subwindows were examined only as sensitivity analyses.

The primary outcome was $Y_i=1$ for a generated delivery interval beginning at gestational week 20 through 36 and $Y_i=0$ when no delivery occurred before week 37. All pregnancies were ongoing at entry and were completely observed for this generated endpoint. The stored value 37 denotes the absence of an earlier event; it is not an exact gestational age for term deliveries and was not analysed as a continuous delivery age.

The generator does not distinguish live births from pregnancy losses or spontaneous from indicated delivery. The endpoint should therefore be read as generated delivery before week 37, rather than a clinically validated preterm-birth diagnosis. The exploratory endpoint used a threshold before week 34 under the same limitations.

### Association models

The prespecified primary model was frequentist logistic regression,

$$
\log\left\{\frac{p_i}{1-p_i}\right\}=\alpha+\beta_H(H_i/5)+\beta_A\{(A_i-30)/5\}+\gamma^{\mathsf T}Z_i,
\qquad p_i=\Pr(Y_i=1\mid H_i,A_i,Z_i).
$$

The adjustment vector contained smoking, deprivation, fixed location indicators, log one plus mean precipitation during the same 14-day window, annual sine and cosine terms at entry, and calendar year centred at 2022. Johannesburg, no smoking and no deprivation were reference categories. Covariates were retained for the specified adjustment structure, rather than selected by stepwise significance testing. Separate univariable models provided the unadjusted estimates. A maternal-adjusted heat model was retained to show the effect of adding location and calendar adjustment.

The exposure odds ratio was $\exp(\beta_H)$ per 5 degrees C; the age odds ratio was per five years. Odds ratios were not interpreted as risk ratios. HC3 sandwich covariance accounted for leverage in the estimated score variability [3]. With $I=X^{\mathsf T}WX$ and logistic leverage $h_i$, the estimator was

$$
\widehat V_{HC3}=I^{-1}\left\{\sum_i\frac{x_ix_i^{\mathsf T}(Y_i-\hat p_i)^2}{(1-h_i)^2}\right\}I^{-1}.
$$

Confidence intervals were $\exp\{\hat\beta\pm1.96\operatorname{SE}_{HC3}(\hat\beta)\}$, with two-sided Wald p values. These are confidence intervals from a frequentist analysis, not posterior credible intervals. Robust standard errors do not remove confounding or correct a misspecified mean relationship.

### Model diagnostics and sensitivity analyses

Diagnostics assessed convergence, finite estimates, model-matrix rank and complete or quasi separation. Separation was checked by linear programming: a nonzero coefficient direction with nonnegative outcome-signed margins and at least one strict margin indicates separation. The implementation was checked on known complete, quasi and overlapping examples.

Heat and age linearity on the logit scale were assessed by separately replacing each linear term with a three-degree-of-freedom natural cubic spline and using a likelihood-ratio test. Normal distributions of predictors or residuals are not assumptions of logistic regression. A spline sensitivity model for both heat and age supplied a pointwise association curve relative to median heat and an exploratory interquartile contrast; this model was not selected because of its p values.

Collinearity was assessed using generalised variance-inflation factors and their degrees-of-freedom adjustment. Influence was examined with leverage and Cook's distance. The primary fit retained every eligible record. A sensitivity analysis omitted only the five largest Cook distances and disclosed the resulting sample and event counts; the screening threshold 4/N was not an exclusion rule.

Grouped goodness of fit was examined using ten fitted-probability groups and a Hosmer-Lemeshow statistic. Its p value is grouping dependent and does not prove that a model is correct. The accompanying plot describes in-sample fit rather than independent prediction performance. Additional sensitivity analyses used the separate seven-day windows and location-specific annual seasonal terms. Secondary and sensitivity findings were not multiplicity-adjusted and were treated as exploratory.

### Tables and visualisation

Descriptive characteristics were summarised as median and interquartile range or count and percentage. The gtsummary package produced the participant and regression tables [4], with a custom tidier retaining HC3 intervals and p values. These objects were converted to native editable Word tables with flextable [5]. The same displayed values were exported for the GitHub README. Figures were created with ggplot2. All statistical calculations and narrative values in this analysis were derived in R.

## Results

### Cohort characteristics and exposure

@@COHORT_RESULTS@@

@@TABLE1@@

Continuous entries are median (25th, 75th percentiles); categorical entries are n (%). Smoking and deprivation rows count Yes. Percentages are descriptive and the balanced location allocation reflects the simulation design. No inferential p values were added to the descriptive table.

The environmental distributions differed across locations, with hotter baseline conditions in Abidjan and substantial seasonal variation elsewhere. These plots describe assigned weather exposures, not geographic patterns of clinical disease.

![Figure 1 Distribution of baseline heat by illustrative weather location](FIGPATH01_exposure_distribution.png)

### Heat and maternal associations

@@ASSOCIATION_RESULTS@@

@@TABLE2@@

Johannesburg is the location reference. Adjustment also includes rainfall, season and calendar year. Location estimates describe this synthetic design, not national clinical risks. The forest plot shows the imprecise heat estimate alongside the positive maternal covariate associations.

![Figure 2 Adjusted odds ratios for the primary exposure and maternal characteristics](FIGPATH02_adjusted_associations.png)

### Functional form and other diagnostics

@@DIAGNOSTIC_RESULTS@@

@@DIAGNOSTICS_TABLE@@

The spline curve is a sensitivity assessment of functional form, with pointwise intervals rather than a simultaneous confidence band. It assumes a common heat response across locations and does not estimate a local curve for every city. The plotted range is restricted to the pooled 1st-99th exposure percentiles; local overlap remains uneven.

@@SPLINE_RESULTS@@

![Figure 3 Exploratory nonlinear heat association relative to the median exposure](FIGPATH03_nonlinear_association.png)

@@INFLUENCE_RESULTS@@

![Figure 4 Influence screening while retaining the full primary cohort](FIGPATH04_influence.png)

The grouped plot provides a descriptive check of observed and fitted event proportions. Agreement is evaluated within the same cohort used to fit the regression and is not a validation of an operational prediction tool.

![Figure 5 Grouped in sample goodness of fit](FIGPATH05_grouped_fit.png)

### Exposure and adjustment sensitivities

@@SENSITIVITY_RESULTS@@

@@SENSITIVITY_TABLE@@

Each heat odds ratio represents a five-degree contrast within the stated window. The before-week-34 analysis uses a different exploratory outcome. The omission analysis removes five events because those records had the five largest Cook distances; it does not replace the full-cohort result. The spline interquartile contrast and linear five-degree contrast are different summaries and are not interchangeable.

## Discussion

The main finding was that the crude positive association between baseline heat and subsequent delivery before week 37 did not persist in the fully adjusted primary linear model. The change in magnitude and direction after location and calendar adjustment shows that the crude estimate was sensitive to the covariate structure. It does not prove that heat is protective or that every source of confounding has been removed.

Maternal age, smoking and deprivation had positive adjusted associations in the synthetic cohort. These relationships are coherent with the specified generator, but their recovery is not new clinical evidence. The environmental and maternal variables also have different interpretations: heat was a short baseline exposure measurement, while the outcome reflected delivery over the following gestational weeks. A baseline association cannot identify the effects of exposure occurring later in pregnancy.

The spline sensitivity suggested that a single linear coefficient does not exhaust the possible exposure pattern. Its interquartile contrast should be interpreted alongside the primary estimate and the exploratory nature of the additional analyses. Isolated intervals from a flexible curve are not sufficient to claim a protective heat effect. The workflow retains these differences rather than selecting whichever result is most statistically significant.

Strengths of the analysis include documented data import and quality checks, a fixed exposure window preceding follow-up, explicit outcome and reference coding, clinically interpretable covariate increments, native gtsummary and flextable outputs, and transparent diagnostics. Missingness and influence were assessed without manufacturing data problems or silently excluding observations. The model was estimable and the diagnostic findings were reported with their limitations.

Important limitations arise from the data themselves. All health records are simulated; clinical definitions, treatment, comorbidity, pregnancy loss and live-birth status are absent. Six balanced illustrative locations are not representative clinical centres. One reanalysis coordinate cannot measure personal temperature, housing, occupational exposure, cooling behaviour or mobility. The short baseline window does not represent all subsequent heat, and seasonal adjustment may remain incomplete. Exposure overlaps with location and season, limiting precision and interpretation of a shared coefficient. The endpoint is an interval-coded generated event, not an independently verified clinical diagnosis.

The study therefore demonstrates a complete applied analysis in R, with substantive association estimates and explicit uncertainty. Application to actual patient data would require a clinical data dictionary, validated obstetric endpoints, appropriate covariate and exposure measurement and a justified confounding strategy. These requirements cannot be met merely by choosing a more complicated statistical model.

## Conclusion

In this synthetic cohort, the adjusted primary analysis did not show clear evidence of an association between heat during the 14 days before week-20 entry and later delivery before week 37. Maternal age, smoking and deprivation showed positive adjusted associations. The crude heat association changed after location and calendar adjustment, and functional-form sensitivity warranted cautious interpretation. The study provides reproducible preprocessing, association estimates, diagnostics and publication-style tables while remaining explicitly a simulation-based applied research exercise.

## Data and code availability

The source weather and synthetic cohort are retained unchanged. The complete cleaned analysis dataset, dictionary, exclusion and missingness audit, regression estimates, model objects, diagnostics and table objects are under results/association. The ordered R scripts and run_association.R reproduce the analysis, figures and manuscript source. Native gtsummary and flextable outputs are retained with the exported tables. Weather attribution remains ERA5 through Open-Meteo under the applicable CC BY terms [2].

The complete previous methods study is frozen under archive/methods-study-v1, including its original README, code, results, figures and Word manuscript. It is historical supplementary work rather than the primary study reported here. No clinical cohort, institutional collaboration or empirical discovery is claimed. The manuscript has not undergone journal peer review.

## References

1. Chersich MF and colleagues. Associations between high temperatures in pregnancy and risk of preterm birth, low birth weight, and stillbirths. BMJ 2020;371:m3811. https://www.bmj.com/content/371/bmj.m3811
2. Open-Meteo. Historical Weather API documentation and data attribution. https://open-meteo.com/en/docs/historical-weather-api
3. Zeileis A. Econometric computing with HC and HAC covariance matrix estimators. Journal of Statistical Software 2004;11(10). https://doi.org/10.18637/jss.v011.i10
4. Sjoberg DD and colleagues. gtsummary regression and summary table documentation. https://www.danieldsjoberg.com/gtsummary/
5. Gohel D. flextable documentation. https://ardata-fr.github.io/flextable-book/
)"
values<-list(
 '@@ABSTRACT_RESULTS@@'=sprintf('**Results.** There were %s events among 6,000 pregnancies (%.1f%%). The crude heat OR per five degrees C was %.2f (95%% CI %.2f-%.2f), compared with adjusted OR %s. Adjusted ORs were %s for age per five years, %s for smoking and %s for deprivation. The model converged without detected separation; sensitivity analyses are reported without replacing the primary fit.',format(sum(d$delivery37),big.mark=','),100*mean(d$delivery37),raw_heat$estimate,raw_heat$conf.low,raw_heat$conf.high,heat,age,smoke,depr),
 '@@COHORT_RESULTS@@'=sprintf('All 6,000 source pregnancies were retained, including 796 delivery events before week 37 (13.3%%) and 411 before week 34 (6.9%%). Entry dates ranged from %s to %s. Median baseline heat was %.1f degrees C (IQR %.1f-%.1f). There were 1,000 pregnancies at each illustrative location. Table 1 describes the cohort by the primary outcome.',prep$entry_min,prep$entry_max,median(d$heat14),quantile(d$heat14,.25),quantile(d$heat14,.75)),
 '@@ASSOCIATION_RESULTS@@'=sprintf('The crude heat OR was %.2f (95%% CI %.2f-%.2f) per five degrees C. After adjustment for maternal characteristics, location, rainfall, season and calendar year, the OR was %s. The adjusted interval included one, providing no clear evidence of a linear heat association in the primary analysis. Adjusted ORs were %s for age per five years, %s for smoking and %s for deprivation. These are differences in odds, not percentage changes in risk.',raw_heat$estimate,raw_heat$conf.low,raw_heat$conf.high,heat,age,smoke,depr),
 '@@DIAGNOSTIC_RESULTS@@'=sprintf('The primary model converged with full rank (%d coefficients), finite estimates and no complete or quasi separation detected by linear programming. The likelihood-ratio p values for nonlinear heat and age terms were %s and %s. These tests did not provide strong evidence against linear logit terms, but do not prove linearity. Heat had a one-degree-of-freedom VIF of %.2f and adjusted GVIF of %.2f; this reflects overlap with the adjustment variables and should not be described as an absence of collinearity. The grouped goodness-of-fit statistic was %.2f on %d degrees of freedom (%s).',report_rank,format_p(nonlinearity$p[1]),format_p(nonlinearity$p[2]),vif$GVIF[vif$term=='heat5'],vif$adjusted_GVIF[vif$term=='heat5'],hl,hl_df,p_text(hl_p)),
 '@@SPLINE_RESULTS@@'=sprintf('The exploratory spline contrast comparing the 75th exposure percentile (%.1f degrees C) with the 25th (%.1f degrees C) gave OR %.2f (95%% CI %.2f-%.2f). This differs from the prespecified constant five-degree linear contrast and is not evidence that higher heat prevents adverse outcomes. Additional contrasts were not corrected for multiplicity.',diag$spline_iqr$high,diag$spline_iqr$low,diag$spline_iqr$OR,diag$spline_iqr$lower,diag$spline_iqr$upper),
 '@@INFLUENCE_RESULTS@@'=sprintf("There were %d records above the Cook distance screening threshold 4/N, with maximum Cook distance %.4f. This low screening threshold was not used to delete records from the primary cohort. Omitting the five largest values yielded OR %.2f (95%% CI %.2f-%.2f), with 5,995 pregnancies and 791 events, compared with the full-cohort OR %.2f.",sum(cooks>4/nrow(d)),max(cooks),influence_row$estimate,influence_row$conf.low,influence_row$conf.high,primary_terms$estimate[primary_terms$term=='heat5']),
 '@@SENSITIVITY_RESULTS@@'=sprintf('The separate seven-day windows retained imprecise fully adjusted estimates. Allowing location-specific seasonal terms gave OR %.2f (95%% CI %.2f-%.2f). For the exploratory before-week-34 endpoint, the OR was %.2f (95%% CI %.2f-%.2f). These checks provide context for the primary association, not independent confirmatory findings.',season_row$estimate,season_row$conf.low,season_row$conf.high,early_row$estimate,early_row$conf.low,early_row$conf.high))
for(n in names(values))body<-gsub(n,values[[n]],body,fixed=TRUE)
table_markers<-c('@@FLOW_TABLE@@'='flow','@@TABLE1@@'='table1','@@TABLE2@@'='table2','@@DIAGNOSTICS_TABLE@@'='diagnostics','@@SENSITIVITY_TABLE@@'='sensitivity')
readme_body<-body;word_body<-body
for(n in names(table_markers)){
 name<-table_markers[[n]];md<-paste(readLines(file.path(tabdir,paste0(name,'.md')),encoding='UTF-8',warn=FALSE),collapse='\n')
 readme_body<-gsub(n,md,readme_body,fixed=TRUE)
 # Return the flextable object itself so knitr dispatches its native DOCX method.
 chunk<-paste0('```{r}\n#| echo: false\n#| message: false\n#| warning: false\nreadRDS("../results/association/tables/',name,'_flextable.rds")\n```')
 word_body<-gsub(n,chunk,word_body,fixed=TRUE)
}
title<-'Maternal heat exposure and subsequent delivery before 37 weeks'
subtitle<-'An association analysis in R using a simulated pregnancy cohort'
repro<-r"(
## Reproducing the R association study

Use R 4.5.1 or a compatible version. Required packages are jsonlite, sandwich, car, lpSolve, gtsummary, flextable, officer, ggplot2, dplyr, broom, knitr and rmarkdown. The installed versions and complete R session are recorded with the association outputs. Quarto is needed only for Word rendering.

```r
install.packages(c('jsonlite','sandwich','car','lpSolve','gtsummary',
                   'flextable','officer','ggplot2','dplyr','broom',
                   'knitr','rmarkdown'), repos='https://cloud.r-project.org')
```

```bash
Rscript --vanilla run_association.R --analysis-only
Rscript --vanilla run_association.R
Rscript --vanilla run_association.R --render
```

The first command runs preprocessing, association models, diagnostics and data-integrity checks. The second also produces gtsummary and flextable objects, ggplot figures and the manuscript source. The third renders the Word manuscript. Set QUARTO_PATH to the Quarto executable if it is not on PATH. The workflow does not regenerate the original cohort or run repeated method comparisons. Rendering uses R with user profile and environment files disabled.

The ordered implementation is R/01_preprocess.R, R/02_associations.R, R/03_tables_figures.R and R/04_manuscript.R. Configuration is in association_config.json. R/validate_association.R checks source preservation, timing, direct weather aggregation, model convergence, exposure-unit scaling and the separation diagnostic. These are software and data checks, not the research objectives.

The original methods study, including its independent R/Python comparisons, is preserved in [archive/methods-study-v1](archive/methods-study-v1/README.md). Its results are not combined with the association findings. The root-level historical scripts remain for provenance; run_association.R is the main entry point.
)"
header<-paste0('# ',title,'\n\n**',subtitle,'**\n\n**Vusumuzi Mabasa | 26 November 2026**\n\n**Public ERA5 weather with entirely simulated maternal characteristics and delivery outcomes.**\n\n[Word manuscript](report/Heat_Health_Manuscript.docx) | [Quarto source](report/Heat_Health_Manuscript.qmd) | [R workflow](run_association.R) | [Association results](results/association) | [Original study archive](archive/methods-study-v1/README.md)\n\n')
writeLines(paste0(header,gsub('FIGPATH','figures/association/',readme_body,fixed=TRUE),repro),file.path(root,'README.md'),useBytes=TRUE)
yaml<-paste0('---\ntitle: "',title,'"\nsubtitle: "',subtitle,'"\nauthor: "Vusumuzi Mabasa"\ndate: "2026-11-26"\ndate-format: "D MMMM YYYY"\nformat:\n  docx:\n    toc: false\n    reference-doc: association_reference.docx\nfontsize: 11pt\nexecute:\n  echo: false\n  warning: false\n  message: false\n---\n\n**Association study with public weather and entirely synthetic health outcomes**\n\n```{r}\nlibrary(flextable)\n```\n\n')
writeLines(paste0(yaml,gsub('FIGPATH','../figures/association/',word_body,fixed=TRUE)),file.path(root,'report/Heat_Health_Manuscript.qmd'),useBytes=TRUE)
writeLines(c('R ASSOCIATION STUDY WORKFLOW','Vusumuzi Mabasa','All health outcomes are simulated.','',repro),file.path(root,'workflow.txt'),useBytes=TRUE)
cat('Research README and Quarto manuscript generated in R\n')
