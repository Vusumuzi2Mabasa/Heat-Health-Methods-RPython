from common import *
import pandas as pd
from scipy.stats import norm
from sklearn.metrics import roc_auc_score
out=ROOT/'results';r=pd.read_csv(out/'r_fits.csv',keep_default_na=False);py=pd.read_csv(out/'python_fits.csv',keep_default_na=False);mf=pd.read_csv(ROOT/'data/simulation_manifest.csv',keep_default_na=False)
assert len(r)==len(py)==1500
assert r.converged.all() and py.converged.all()
keys=['replicate','scenario','model'];both=r.merge(py,on=keys,suffixes=('_r','_python'),validate='one_to_one')
differences={k:float((both[k+'_r']-both[k+'_python']).abs().max()) for k in ['estimate','se','loglik']}
assert differences['estimate']<1e-5 and differences['se']<1e-4 and differences['loglik']<1e-5,differences
co=pd.read_csv(out/'r_fits_coefficients.csv',keep_default_na=False).merge(pd.read_csv(out/'python_fits_coefficients.csv',keep_default_na=False),on=['scenario','model','term'],suffixes=('_r','_python'),validate='one_to_one')
differences['all_coefficients_replicate1']=float((co.estimate_r-co.estimate_python).abs().max())
assert differences['all_coefficients_replicate1']<1e-4
rows=[]
for (scenario,model),g in r.groupby(['scenario','model']):
 cover=(g.lower<=g.truth)&(g.upper>=g.truth);reject=g.p<.05;R=len(g);p=reject.mean();den=1+1.96**2/R;center=(p+1.96**2/(2*R))/den;half=1.96*np.sqrt(p*(1-p)/R+1.96**2/(4*R**2))/den
 rows.append(dict(scenario=scenario,model=model,replicates=R,truth=g.truth.iloc[0],mean_estimate=g.estimate.mean(),bias=(g.estimate-g.truth).mean(),bias_mcse=g.estimate.std()/np.sqrt(R),empirical_sd=g.estimate.std(),mean_se=g.se.mean(),rmse=np.sqrt(np.mean((g.estimate-g.truth)**2)),coverage=cover.mean(),coverage_mcse=np.sqrt(cover.mean()*(1-cover.mean())/R),rejection_rate=p,rejection_mcse=np.sqrt(p*(1-p)/R),rejection_lower=center-half,rejection_upper=center+half,mean_events=g.events.mean(),mean_observed_events=g.observed_events.mean()))
pd.DataFrame(rows).to_csv(out/'simulation_summary.csv',index=False)
added=[];rng=np.random.default_rng(CONFIG['seed']+333)
for language,path in [('R','r_landmark_predictions.csv'),('Python','landmark_predictions.csv')]:
 pred=pd.read_csv(out/path)
 for test,g in pred.groupby('test'):
  baseline=g[g.model=='Maternal baseline'][['id','outcome','prediction']].rename(columns={'prediction':'baseline'})
  for model in ['Logistic','Random forest','Gradient boosting']:
   d=baseline.merge(g[g.model==model][['id','prediction']],on='id',validate='one_to_one');y=d.outcome.to_numpy();a=d.baseline.to_numpy();b=d.prediction.to_numpy();delta=[]
   for k in range(300):
    j=rng.integers(0,len(d),len(d))
    if len(np.unique(y[j]))==2:delta.append(roc_auc_score(y[j],b[j])-roc_auc_score(y[j],a[j]))
   lo,hi=np.quantile(delta,[.025,.975]);added.append(dict(language=language,test=test,model=model,delta_auc=roc_auc_score(y,b)-roc_auc_score(y,a),lower=lo,upper=hi,delta_brier=np.mean((y-b)**2)-np.mean((y-a)**2)))
pd.DataFrame(added).to_csv(out/'added_value.csv',index=False)
# Independent logistic coefficient cross-check.
lr=pd.read_csv(out/'r_landmark_logistic.csv');lp=pd.read_csv(out/'python_landmark_logistic.csv');diff=float((lr.coefficient-lp.coefficient).abs().max());assert diff<1e-5
checks=dict(fit_records_per_language=len(r),converged_r=int(r.converged.sum()),converged_python=int(py.converged.sum()),max_absolute_differences=differences,landmark_logistic_max_coefficient_difference=diff,weather_rows=13152,independent_cohorts=100,scenario_datasets=500,oracle_reuse='Same fit as nonlinear in null, moderate and nonlinear scenarios because exposure and rows are identical',scope='Numerical reproducibility validates implementation, not clinical validity')
(out/'validation.json').write_text(json.dumps(checks,indent=2))
pd.DataFrame([dict(scenario=s,mean_events=g.events.mean(),min_events=g.events.min(),max_events=g.events.max(),mean_event_fraction=(g.events/g.n).mean()) for s,g in mf.groupby('scenario')]).to_csv(out/'cohort_summary.csv',index=False)
print(json.dumps(checks,indent=2))
