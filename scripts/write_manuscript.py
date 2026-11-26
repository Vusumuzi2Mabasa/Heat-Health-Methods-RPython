from pathlib import Path
import json
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
s=pd.read_csv(ROOT/'results/simulation_summary.csv',keep_default_na=False);m=pd.read_csv(ROOT/'results/r_prediction_metrics.csv');mp=pd.read_csv(ROOT/'results/prediction_metrics.csv');a=pd.read_csv(ROOT/'results/added_value.csv');splits=pd.read_csv(ROOT/'results/prediction_splits.csv');v=json.loads((ROOT/'results/validation.json').read_text());std=json.loads((ROOT/'results/standardisation.json').read_text())
support=pd.read_csv(ROOT/'results/exposure_support.csv')
def row(sc,mo='nonlinear'):return s[(s.scenario==sc)&(s.model==mo)].iloc[0]
def table(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])
def label(x):return {'null':'Null','moderate':'Moderate','nonlinear':'Nonlinear','exposure_error':'Exposure error','informative_missingness':'Informative missingness','linear':'Current linear','oracle':'Oracle','temporal':'Temporal','geographic':'Geographic','combined_shift':'Combined shift','train':'Training','excluded_transition':'Excluded transition'}.get(x,x)
def figure(name,caption):return f'![{caption}](FIGURES/{name}.png)\n'
title='Heat exposure and preterm birth across African settings'
subtitle='A simulation study of delayed effects exposure uncertainty and predictive transportability'
null=row('null');moderate=row('moderate');nonlinear=row('nonlinear');err=row('exposure_error');missing=row('informative_missingness')
simtable=table(['Scenario','Model','Bias','Coverage','Rejection','MCSE'],[[label(r.scenario),label(r.model),f'{r.bias:+.3f}',f'{100*r.coverage:.0f}%',f'{100*r.rejection_rate:.0f}%',f'{100*r.rejection_mcse:.1f}%'] for r in s.itertuples()])
predtable=table(['Test set','R model','AUC 95% CI','Brier','Slope'],[[label(r.test),r.model,f'{r.auc:.3f} ({r.auc_lower:.3f}–{r.auc_upper:.3f})',f'{r.brier:.3f}',f'{r.calibration_slope:.2f}'] for r in m.itertuples()])
splittable=table(['Partition','Pregnancies','Preterm births'],[[label(r.split),r.n,r.events] for r in splits.itertuples()])
addedtable=table(['Test set','Model','Change in AUC 95% CI','Change in Brier'],[[label(r.test),r.model,f'{r.delta_auc:+.3f} ({r.lower:+.3f} to {r.upper:+.3f})',f'{r.delta_brier:+.4f}'] for r in a[a.language=='R'].itertuples()])
supporttable=table(['Location','Current window range °C','Both windows above 25°C','Both windows at least 30°C'],[[r.site.replace('_',' '),f'{r.min_T0:.1f}–{r.max_T0:.1f}',f'{r.fraction_both_above25:.1%}',f'{r.fraction_both_atleast30:.1%}'] for r in support.itertuples()])
body=r'''## Abstract

**Background.** Heat and pregnancy research requires methods that respect gestational risk sets, delayed exposure effects and the distinction between association and prediction. This study evaluated such methods using public historical weather and entirely synthetic pregnancies.

**Methods.** ERA5 reanalysis weather from six African locations during 2019–2024 informed 100 independently seeded cohorts of 1,800 pregnancies. Five paired scenarios represented a null effect, a moderate heat effect, nonlinear delayed effects, unobserved personal exposure and outcome-dependent missingness. Complementary log-log weekly survival models were implemented independently in R and Python. An additional 6,000-pregnancy cohort supported week-28 prediction with maternal-only ridge logistic regression, weather-augmented ridge regression, Random Forest and gradient boosting. Temporal, geographic and combined test sets were held out.

ABSTRACT_RESULTS

**Conclusions.** This is an implementation and performance study under known simulated mechanisms. It does not estimate clinical effects in African populations or validate an operational heat-warning system. Performance depends on exposure measurement, missingness and the prediction setting; additional model complexity does not guarantee better prediction.

## Introduction

Ambient heat is a relevant maternal and newborn health exposure. A systematic review reported associations between higher temperatures and adverse pregnancy outcomes, including preterm birth [1]. These findings motivate careful study designs; they do not establish the numerical parameters used in the present simulations. Wits Planetary Health Research describes programmes addressing heat, maternal health and early warning, including HEAT Center and Bio-HEAT [2,3]. This independent methodological portfolio study is motivated by those topics and has no affiliation with, access to, or validation against their participant data.

Three statistical questions are distinguished. First, can a time-indexed model recover a known heat association without using exposure after delivery? Second, how does recovery change when ambient weather differs from personal exposure or when exposure records are selectively missing? Third, can information available at a pregnancy landmark predict later preterm delivery in different years and locations? A final explanatory exercise examines model-standardised risks under a hypothetical reduction in temperature. The study was designed to evaluate these questions rather than maximise the apparent performance of one algorithm.

## Methods

### Design and data provenance

Daily maximum and minimum air temperature and precipitation were obtained through the Open-Meteo Historical Weather API with explicit ERA5 model selection [4]. The locations were Johannesburg, Durban, Cape Town, Harare, Nairobi and Abidjan. Each location contributed 2,192 consecutive days from 1 January 2019 to 31 December 2024, giving 13,152 site-days. API queries, source hashes, returned coordinates and time zones are retained in the provenance directory. The source is gridded reanalysis; it does not measure individual heat exposure. These illustrative city coordinates are not recruitment sites.

Checks required complete consecutive dates, finite temperatures and precipitation, minimum temperature no greater than maximum temperature, plausible temperature bounds and nonnegative rainfall. An initial ERA5-Land query returned unavailable precipitation; the final study therefore uses a consistent ERA5 product across all sites, without filling the missing weather values. Monthly weather summaries describe the exposure environment but are not evidence of changes in health outcomes.

WEATHER_FIGURE

All participants, covariates and birth outcomes were generated synthetically. The original heat–pregnancy project was preserved in its separate directory. The present design does not reproduce patient data, recruitment or effect estimates from Wits PHR. No clinical records, names or participant identifiers are included. Original code and manuscript are authored by Vusumuzi Mabasa.

### Pregnancy risk sets and temporal ordering

Each simulated cohort contained 1,800 pregnancies balanced across the six locations. Entry at the beginning of gestational week 20 was sampled between 15 January 2020 and 1 July 2024. Age followed a normal distribution with mean 29 and standard deviation 6 years, clipped to 18–44. Smoking and deprivation were independent Bernoulli variables with probabilities 0.15 and 0.40. Independence is a simplifying assumption rather than a claim about observed social patterns.

Follow-up comprised weekly intervals beginning at gestational weeks 20 through 36. A pregnancy contributed an interval only while undelivered; follow-up stopped at the first simulated delivery. Delivery by the end of week 36 was the event. Pregnancies with no earlier event were assigned week 37 for this simulation. Miscarriage, stillbirth, spontaneous versus indicated delivery and loss to follow-up were not modelled. The target population is therefore hypothetical pregnancies still ongoing at week 20, not all conceptions.

At interval entry on day $d$, two exposure windows used only previously completed days:

$$
T_{i0}(d)=\frac{1}{7}\sum_{j=1}^{7}T_{s_i}(d-j),\qquad
T_{i1}(d)=\frac{1}{7}\sum_{j=8}^{14}T_{s_i}(d-j).
$$

Here $T_s$ is daily maximum temperature at location $s$. Mean precipitation over days 1–7, calendar-year trend and annual sine and cosine terms were included as covariates. Person-week records after delivery were excluded. These rules avoid using a whole-pregnancy exposure average whose duration is determined by delivery time.

### Generating model and estimand

Let $h_{ig}$ denote the probability of delivery in week $g$, conditional on remaining undelivered. The generating model was

$$
\log[-\log(1-h_{ig})]=\alpha_{s_i}+0.84u_g+\gamma^\top X_{ig}
+\sum_{\ell=0}^{1}\{\beta_{\ell 1}z_{i\ell g}+\beta_{\ell 2}z_{i\ell g}^{2}\}.
$$

$$
u_g=(g-28)/4,\qquad z_{i\ell g}=\max(T_{i\ell g}-25,0)/5.
$$

The reference intercept was −5.62, with five fixed site offsets 0.10, −0.10, 0.05, −0.05 and 0.12. Covariate coefficients were 0.12 for age centred at 30 and divided by 10, 0.40 for smoking, 0.25 for deprivation, −0.025 for log-transformed mean rainfall, 0.15 and −0.10 for seasonal sine and cosine, and 0.025 per calendar year centred at 2022. A quadratic gestational-age term was estimated with generating coefficient zero. Parameters were prespecified illustrative assumptions; neither prevalence nor heat effects were estimated from a patient cohort.

The nonlinear heat coefficients were $(0.18,0.08,0.12,0.04)$. They were halved in the moderate scenario and set to zero in the null scenario. The primary estimand compared 30°C with 25°C in **both** completed exposure windows:

$$
\theta=\beta_{01}+\beta_{02}+\beta_{11}+\beta_{12},\qquad
H=\exp(\theta).
$$

True $\theta$ was 0, 0.21 or 0.42, corresponding to $H=1$, 1.234 or 1.522. Under the complementary log-log model, $H$ is the ratio of underlying interval cumulative hazards conditional on covariates. It is not a cumulative preterm-birth risk ratio. The fixed 25°C threshold is an illustrative basis knot, not a demonstrated clinical threshold or a location-specific adaptation threshold.

### Scenarios and model comparisons

One hundred base cohorts were generated with independent seeds; the same underlying draws were reused across the five scenarios within each replicate. Consequently, there are 500 scenario datasets but 100 independent base replicates. All scenarios share a fixed historical weather panel. The exposure-error scenario added a person-specific $U_i\sim N(0,2^2)$ degrees to both ambient windows when generating outcomes. The analyst observed ambient weather only; the oracle knew $U_i$. This is a stable personal-versus-ambient offset, resembling Berkson-type exposure error, rather than a fully realistic wearable measurement process.

In the informative-missingness scenario the probability that interval exposure was unavailable was

$$
\operatorname{logit}\Pr(M_{ig}=1)=-3.2+0.45z_{i0g}+0.5D_i+1.5Y_{ig}.
$$

Here $D_i$ denotes deprivation and $Y_{ig}$ indicates delivery in that interval. Birth outcomes continued to define follow-up even if an exposure row was dropped. This intentionally outcome-dependent mechanism is not missing at random. No imputation or weighting correction is claimed.

The observed nonlinear model estimated the four heat terms and 15 baseline terms, including the intercept. A deliberately misspecified comparator retained the same baseline terms but used current-window temperature linearly, omitting the delayed window and nonlinear basis. Its single coefficient was evaluated against the joint heat estimand as an illustration of the consequences of using an inadequate model; it does not estimate a distinct correctly specified target. The oracle used true personal exposures and complete risk sets. Oracle fits equal the nonlinear fits in the null, moderate and nonlinear scenarios and were reused there rather than refitted.

This finite exposure-by-window basis follows the general principle of distributed lag nonlinear modelling [5], but it is a custom hinge-quadratic basis with two lag-window indicators. It is not a fitted spline, a full daily lag surface or a run of the R dlnm package. Seasonal and site adjustment reproduce the specified data-generating structure; they do not demonstrate control of all real-world confounding.

### Estimation and simulation performance

R used binomial complementary log-log GLMs; Python used statsmodels GLMs. Each language independently reconstructed exposure windows, risk sets and model matrices from the same cohort-level files. Maximum likelihood used a maximum of 100 iterations and tolerance $10^{-10}$. Coefficient covariance matrices yielded

$$
\widehat{\operatorname{se}}(\hat\theta)=\sqrt{c^\top\widehat Vc},\qquad
CI_{95}=\hat\theta\pm1.96\widehat{\operatorname{se}}(\hat\theta).
$$

The vector $c$ selects the four heat coefficients. A two-sided Wald test of $\theta=0$ used a nominal 5% level. This is one prespecified contrast, not an omnibus four-parameter heat test. Bias, empirical standard deviation, mean standard error, root mean squared error, coverage and rejection rate were recorded. Rejection is type I error in the null scenario and power in the other scenarios. For a proportion $\hat p$ over $R=100$ repetitions,

$$
\operatorname{MCSE}(\hat p)=\sqrt{\hat p(1-\hat p)/R},\qquad
\operatorname{MCSE}(\widehat{\mathrm{bias}})=s_{\hat\theta}/\sqrt R.
$$

The maximum possible proportion MCSE is 0.05. Wilson intervals for rejection rates and all replicate estimates are supplied as results files. The design supports a manageable methodological demonstration, not highly precise power calculations. Model-based standard errors are appropriate under the correctly specified generator; their behaviour under misspecification is itself part of the evaluation.

### Landmark prediction and independent testing

An independently seeded illustrative cohort contained 6,000 pregnancies. Among those still undelivered at the beginning of week 28, the prediction endpoint was subsequent delivery before week 37. Predictors comprised age, smoking, deprivation, two prior temperature windows, recent rainfall, annual seasonal terms and calendar year at the landmark. No future weather, latent personal offset or post-delivery information was available to the classifiers.

Training used Johannesburg, Durban, Cape Town and Harare, restricted to pregnancies whose entire planned 119-day follow-up from week-20 entry ended by 31 December 2022. The temporal test used entries on or after 1 January 2023 in those locations. Nairobi and Abidjan were completely withheld from training: their earlier entries formed the geographic test and their later entries the combined-shift test. Pregnancies in the intervening calendar gap were excluded from prediction evaluation. Geographic validation here concerns two illustrative locations, not general validation across Africa.

PREDICTION_SPLITS

The maternal baseline used ridge logistic regression with three maternal covariates. An augmented ridge model used all nine predictors. With standardisation estimated in training, ridge minimises

$$
-\sum_{i\in\mathcal T}\{y_i\log p_i+(1-y_i)\log(1-p_i)\}+\tfrac12\|b\|_2^2,
\qquad p_i=\operatorname{logit}^{-1}(b_0+x_i^\top b).
$$

The intercept was unpenalised. Python used $C=1$; R used glmnet with $\alpha=0$ and $\lambda=1/n_{\mathcal T}$ to match the loss normalisation. Random Forest averaged probabilities from 400 bootstrap trees, with minimum terminal size 25 and six candidate predictors per split. Gradient boosting added 150 depth-2 trees at learning rate 0.03, minimum node size 25 and Bernoulli deviance:

$$
F_m(x)=F_{m-1}(x)+0.03f_m(x),\qquad
p(x)=\operatorname{logit}^{-1}\{F_{150}(x)\}.
$$

The trees approximate successive negative gradients of the binomial loss. Settings were prespecified; no test-set tuning was performed. R ranger and gbm and Python scikit-learn implement different tree algorithms, so identical tree predictions are neither required nor asserted. R prediction is primary in the manuscript, with complete Python results retained. An additional unpenalised logistic specification was cross-checked numerically between languages.

Evaluation included ROC AUC, average precision, Brier score, calibration intercept with logit predictions as an offset and calibration slope from a logistic regression on those logits. Ideal intercept and slope are 0 and 1. Participant bootstrap percentile intervals used 300 resamples, conditional on the fitted model and test sample structure. Paired resampling estimated changes in AUC against the maternal baseline. These intervals do not include training-set uncertainty or uncertainty across future cohorts. Binned calibration plots are descriptive and are not formal calibration tests.

### Hypothetical exposure standardisation

For the illustrative cohort, fitted weekly hazards were combined into cumulative risk,

$$
R_i(a)=1-\prod_{g=20}^{36}\{1-\hat h_{ig}(T_{i0g}-a,T_{i1g}-a)\},\qquad
\Delta=\frac1N\sum_i\{R_i(2)-R_i(0)\}.
$$

The exercise compared ambient exposure with a hypothetical uniform 2°C reduction. Full latent weather schedules were used only for this retrospective explanatory calculation and for its known generating truth, never as prediction inputs. An analytic delta-method interval used the gradient of the mean risk difference with respect to the fitted coefficients, with variance $\nabla\Delta^\top\widehat V\nabla\Delta$. It excludes uncertainty in the population, weather and intervention mechanism. This calculation is not an evaluation of a feasible cooling intervention, a causal effect identified in patient data or attribution to anthropogenic climate change.

## Results

### Simulation performance and numerical reproducibility

SIMULATION_NARRATIVE

SIMULATION_TABLE

Bias refers to the log joint heat contrast; coverage is nominally 95%. The rejection MCSE is a percentage-point uncertainty measure across 100 replicates. Complete results, including coverage MCSE and Wilson intervals, are provided in results/simulation_summary.csv. The linear comparator is deliberately misspecified and oracle results are an upper-information benchmark.

The bias comparison separates the error introduced by an insufficient exposure model from the error introduced by unavailable information. The Monte Carlo intervals describe uncertainty in average bias; they are not confidence intervals for individual pregnancies.

BIAS_FIGURE

The calibration plot should be interpreted jointly with its finite simulation precision. A rejection rate close to 5% in a null scenario is reassuring for this generator, but 100 replicates cannot establish precise tail calibration for all plausible heat-health datasets.

SIMULATION_CALIBRATION_FIGURE

### Illustrative exposure response

The illustrative curve compares equal temperatures in both lag windows with 25°C in both windows. Its uncertainty reflects a four-term contrast and is pointwise, without simultaneous-band correction. Agreement with the planted curve evaluates recovery in this synthetic cohort; divergence at hot values can reflect limited exposure support and coefficient uncertainty. A flat fitted curve below 25°C is imposed by the chosen hinge basis and is not an empirical discovery of no cold effect.

The 30°C joint exposure occurred predominantly in Abidjan. Some locations had no risk intervals with both windows at least 30°C, while Abidjan had no current-window mean as low as 25°C. The pooled contrast therefore relies on the assumed common response across sites and should not be interpreted as a separately identified local heat effect in every city. The curve extends to 35°C, slightly beyond the observed current-window range, to show the specified basis; its upper edge is extrapolation. The association curve and standardisation model use the full illustrative cohort; their explanatory fit is separate from the held-out prediction exercise.

EXPOSURE_SUPPORT_TABLE

RESPONSE_FIGURE

### Held out preterm prediction

PREDICTION_NARRATIVE

PREDICTION_TABLE

The calibration slope measures the spread of predictions relative to observed outcomes: slopes below 1 suggest predictions that vary too strongly, while slopes above 1 suggest insufficient variation, conditional on this evaluation. These are held-out synthetic results; none constitutes external clinical validation.

ROC_FIGURE

CALIBRATION_FIGURE

ADDED_VALUE_NARRATIVE

ADDED_VALUE_TABLE

Positive changes in AUC favour the augmented model; negative changes in Brier score favour it. The augmented models add weather, seasonal and calendar information together, so these comparisons do not isolate temperature alone. Selecting the largest observed AUC after viewing the test results would introduce additional selection bias; no universally preferred algorithm is claimed.

### Hypothetical exposure reduction

STANDARDISATION_NARRATIVE

## Discussion

The principal contribution is a transparent workflow connecting exposure construction, gestational survival analysis, simulation calibration and honest prediction evaluation. The survival and prediction analyses answer different questions. The former estimates a prespecified conditional heat contrast under an explicit model; the latter assesses whether information available at week 28 helps distinguish later events. Useful association recovery need not imply high predictive discrimination, particularly where heat explains only a small part of outcome variation and future exposure remains unknown.

Low power in this design should not be interpreted as evidence that heat has no effect: the true effect is known because it was generated. Each cohort provides far fewer events than participants, the two exposure windows are correlated, and hot exposure support is uneven across locations. These features limit precision even when the fitted model is correct. Changing the simulated effect or selecting a favourable replicate to obtain a significant result would answer a different question. The reported configuration and results are retained as generated.

The oracle is valuable because it distinguishes avoidable modelling limitations from missing information. It knows true personal exposures and retains complete exposure records. Its advantage in a stress scenario does not mean an applied analyst can attain the same result merely by choosing a more flexible model. Exposure measurement, missingness mechanisms and reliable clinical follow-up are substantive requirements. Likewise, a more complex classifier can perform worse in new years or locations even when its training fit is stronger.

Several features improve relevance to environmental epidemiology: real historical weather preserves seasonality and temporal dependence; exposure windows end before each risk interval; pregnancies leave the risk set at delivery; the prediction calendar separates training outcome availability from future evaluation; and geographic testing withholds entire locations. Numerical R–Python agreement provides evidence of consistent implementation. It does not validate the assumed biology, exposure threshold or effect sizes.

Important limitations remain. Each city is represented by one gridded location and no within-city housing, workplace, mobility, humidity, indoor heat or behavioural adaptation process. Covariates are simplified and independent; there is no measured clinical severity or clinical care process. The generator shares its nonlinear basis with one fitted model, favouring that model under correct specification. The misspecified linear comparator is a sensitivity benchmark rather than a claim of superiority over all distributed-lag methods. The homogeneous heat effect, fixed 25°C threshold and two-window memory omit realistic acclimatisation and longer gestational susceptibility. Correlated outcomes due to shared unmeasured site-day causes, residual confounding, competing pregnancy outcomes and informative entry are absent. Complete-case exposure analysis is intentionally inadequate under outcome-dependent missingness.

Only 100 independent base repetitions were conducted and prediction was evaluated in one synthetic cohort. Monte Carlo uncertainty and conditional bootstrap intervals must therefore accompany performance rankings. The standardisation exercise assumes exposure can be shifted without changing other determinants and uses complete simulated schedules; it cannot establish clinical preventability. Personal offsets and missingness parameters are assumptions, not estimated features of a particular cohort.

With appropriately governed patient data, the next steps would be to define the target population and pregnancy time origin, audit delivery and censoring dates, assess exposure coverage and address spatial resolution. A prespecified analysis could compare flexible temperature and lag bases, investigate adaptation using location-specific temperature distributions, model clinically relevant competing events and evaluate exposure error using personal measurements where available. Prediction requires an operational forecast horizon, information available at deployment and prospective recalibration. The current study demonstrates methods that could support such work while making no claim that these extensions have already been validated.

## Conclusion

Public historical weather and synthetic pregnancies provide a reproducible environment for assessing heat-health analysis methods. A correctly specified lagged nonlinear model, an oracle benchmark and prespecified misspecification scenarios make performance interpretable against known truth. Held-out prediction and hypothetical standardisation complement association analysis but do not replace independent clinical validation. The resulting R and Python workflow is a methodological portfolio study, not evidence of patient-level effects or readiness for a public-health warning system.

## Data and code availability

The compact weather panel, synthetic illustrative cohort, configuration, simulation manifest, coefficient estimates, replicate performance, held-out predictions, figures and report source are included. Large raw API responses and regenerable cohort caches are excluded from version control. Source URLs and hashes are retained. Weather attribution: ERA5 reanalysis obtained through Open-Meteo, with source data retaining their applicable CC BY 4.0 attribution requirements [4]. Original code and manuscript use the MIT licence. No Wits PHR participant data were used. No institutional affiliation, endorsement, ethics approval or peer review is implied.

Python scientific dependencies are pinned, R session information and package versions are retained, and a smoke test reproduces matched R and Python fits. The full pipeline regenerates the simulation results before generating figures and the Quarto manuscript. The Word report contains editable native equations rendered by Pandoc. Publication or submission should follow independent scientific review and any necessary institutional approvals.

## References

1. Chersich MF and colleagues. Associations between high temperatures in pregnancy and risk of preterm birth, low birth weight, and stillbirths. BMJ 2020;371:m3811. https://www.bmj.com/content/371/bmj.m3811
2. Wits Planetary Health Research. HEAT Center. https://witsphr.org/Projects/HEAT-Center
3. Wits Planetary Health Research. Bio-HEAT. https://witsphr.org/Projects/Bio-HEAT
4. Open-Meteo. Historical Weather API documentation and data attribution. https://open-meteo.com/en/docs/historical-weather-api
5. Gasparrini A, Armstrong B, Kenward MG. Distributed lag non-linear models. Statistics in Medicine 2010;29:2224–2234. https://pmc.ncbi.nlm.nih.gov/articles/PMC2998707/
'''
abstract_results=f"**Results.** All {v['fit_records_per_language']:,} fit records per language converged. The correctly specified nonlinear model had {null.rejection_rate:.0%} null rejection and {nonlinear.coverage:.0%} interval coverage in the nonlinear scenario. R model AUCs ranged from {m.auc.min():.3f} to {m.auc.max():.3f} across held-out tests; this represents modest discrimination in the simulated prediction problem."
narrative=f"All 1,500 fit records in each language converged. The maximum absolute R–Python differences were {v['max_absolute_differences']['estimate']:.2g} for the joint heat estimate and {v['max_absolute_differences']['se']:.2g} for its standard error. The additional unpenalised landmark logistic model agreed to {v['landmark_logistic_max_coefficient_difference']:.2g} in its coefficients. Oracle records in three scenarios reuse identical nonlinear fits, so these totals do not represent 1,500 distinct numerical optimisations per language.\n\nFor the observed nonlinear model, null rejection was {null.rejection_rate:.0%} (MCSE {null.rejection_mcse:.1%}). Power was {moderate.rejection_rate:.0%} under the moderate contrast and {nonlinear.rejection_rate:.0%} under the full nonlinear contrast. Coverage in those two scenarios was {moderate.coverage:.0%} and {nonlinear.coverage:.0%}. Bias was {err.bias:+.3f} under unobserved personal exposure and {missing.bias:+.3f} under informative exposure missingness. These findings concern the prespecified joint contrast and simulated assumptions, not an estimated clinical heat effect."
rf=m[(m.model=='Random forest')&(m.test=='combined_shift')].iloc[0];mb=m[(m.model=='Maternal baseline')&(m.test=='combined_shift')].iloc[0]
pn=f"The illustrative cohort contained {std['preterm_events']:,} preterm events among 6,000 pregnancies ({std['preterm_events']/6000:.1%}), of which {std['eligible']:,} were still undelivered at week 28. Across the held-out partitions, R AUCs ranged from {m.auc.min():.3f} to {m.auc.max():.3f}; corresponding Python AUCs ranged from {mp.auc.min():.3f} to {mp.auc.max():.3f}. These values indicate modest separation between later preterm and non-preterm outcomes. Tree models did not provide an established general advantage.\n\nFor example, the R Random Forest in the combined-shift test had AUC {rf.auc:.3f}, calibration slope {rf.calibration_slope:.2f} and Brier score {rf.brier:.3f}, compared with maternal-baseline Brier score {mb.brier:.3f}. Ranking ability and probability accuracy should therefore be assessed separately. The test-set sizes and event counts should be considered alongside AUC intervals, Brier scores and calibration slopes."
ar=a[a.language=='R'];an=f"Compared with maternal-only ridge regression, the changes in R AUC ranged from {ar.delta_auc.min():+.3f} to {ar.delta_auc.max():+.3f}. Paired bootstrap intervals are shown below rather than treating every positive change as evidence of benefit. Added predictive value is specific to this synthetic cohort and may differ with another training sample or data-generating mechanism."
sn=f"The fitted model-standardised average risk was {std['baseline_risk']:.2%} under ambient schedules and {std['cooling_risk']:.2%} under the assumed 2°C reduction. The estimated difference was {100*std['risk_difference']:+.2f} percentage points (approximate parameter-only 95% interval {100*std['lower']:+.2f} to {100*std['upper']:+.2f}). The known generating difference was {100*std['truth_risk_difference']:+.2f} percentage points. These values describe the assumed synthetic mechanism and a uniform mathematical exposure shift; they are not estimates of benefits from any real intervention."
replacements={'ABSTRACT_RESULTS':abstract_results,'WEATHER_FIGURE':figure('weather','Figure 1 Public historical weather used to generate exposure histories'),'PREDICTION_SPLITS':splittable,'SIMULATION_NARRATIVE':narrative,'SIMULATION_TABLE':simtable,'BIAS_FIGURE':figure('simulation_bias','Figure 2 Bias in the simulated joint heat contrast'),'SIMULATION_CALIBRATION_FIGURE':figure('simulation_calibration','Figure 3 Coverage and rejection across repeated simulations'),'EXPOSURE_SUPPORT_TABLE':supporttable,'RESPONSE_FIGURE':figure('response','Figure 4 Estimated and generating heat response in the illustrative synthetic cohort'),'PREDICTION_NARRATIVE':pn,'PREDICTION_TABLE':predtable,'ROC_FIGURE':figure('roc','Figure 5 Held out ROC curves from R models'),'CALIBRATION_FIGURE':figure('prediction_calibration','Figure 6 Held out calibration of R models'),'ADDED_VALUE_NARRATIVE':an,'ADDED_VALUE_TABLE':addedtable,'STANDARDISATION_NARRATIVE':sn}
for k,val in replacements.items():body=body.replace(k,val)
repro=r'''
## Reproducing the study

Use Python 3.12, R 4.5.1 and Quarto. R packages used for figures and prediction are ggplot2, ranger, gbm, glmnet and pROC; their installed versions are recorded in provenance. The main survival analysis uses base R. The included ERA5 CSV allows network-free analysis; `--download-weather` refreshes the source through the public API and re-runs download QC.

```bash
python -m venv .venv
# Linux or macOS
source .venv/bin/activate
# Windows PowerShell instead: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
Rscript --vanilla -e "install.packages(c('ggplot2','ranger','gbm','glmnet','pROC'), repos='https://cloud.r-project.org')"
python run_workflow.py --smoke
python run_workflow.py --render
```

On Windows, supply `--rscript "C:/Program Files/R/R-4.5.1/bin/Rscript.exe"` and `--quarto "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe"` if these executables are not on PATH. Full regeneration is substantially slower than the smoke test. Keep R startup files disabled with `--vanilla`; the workflow also disables user environment/profile files for report rendering. Do not overwrite prior releases when changing scenarios: copy configuration and results into a versioned release directory first.

The full pipeline runs `generate.py`, independent `analyse.py` and `analyse.R`, the prediction scripts, invariant tests, numerical comparisons, `figures.R`, and `write_manuscript.py`. These scripts include the full implementation; the report is not a substitute for the executable code. Continuous integration is configured for a two-replicate smoke test. A configured CI file does not imply a remote CI run has been observed.

Outputs: `results/simulation_summary.csv` contains all simulation metrics; `results/r_prediction_metrics.csv` and `results/prediction_metrics.csv` contain R and Python prediction results; `results/added_value.csv` contains paired AUC comparisons; `results/validation.json` records numerical checks; `report/Heat_Health_Manuscript.docx` is the Word manuscript. All health data and reported health effects are simulated.
'''
(ROOT/'README.md').write_text(f'# {title}\n\n**{subtitle}**\n\n**Vusumuzi Mabasa · 26 November 2026**\n\n**Simulation study. Public weather with entirely synthetic pregnancies and outcomes.**\n\n'+body.replace('FIGURES/','figures/')+repro,encoding='utf-8')
qmd=f'''---
title: "{title}"
subtitle: "{subtitle}"
author: "Vusumuzi Mabasa"
date: "26 November 2026"
date-format: "D MMMM YYYY"
format:
  docx:
    toc: false
    number-sections: false
fontsize: 11pt
execute:
  echo: false
  warning: false
  message: false
---

**Simulation study using public weather and entirely synthetic pregnancies and outcomes**

```{{r}}
stopifnot(file.exists('../results/validation.json'))
```

'''+body.replace('FIGURES/','../figures/')
report=ROOT/'report';report.mkdir(exist_ok=True);(report/'Heat_Health_Manuscript.qmd').write_text(qmd,encoding='utf-8')
(ROOT/'workflow.txt').write_text('HEAT HEALTH METHODS R AND PYTHON\nAuthor Vusumuzi Mabasa\nSIMULATED HEALTH DATA ONLY\n\n'+repro.replace('```bash','').replace('```','')+'\n\n'+(ROOT/'results/validation.json').read_text(),encoding='utf-8')
print('Saved research manuscript source and README')
