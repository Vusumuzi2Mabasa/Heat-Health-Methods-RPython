from common import *
import pandas as pd
import statsmodels.api as sm
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier,GradientBoostingClassifier
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss,roc_curve
from scipy.special import logit
import importlib.metadata as metadata
w=weather();c=load_cohort(ROOT/'data/illustrative_cohort.csv');h=histories(c,w);n=len(c['id']);p=h['person'];birth=c['birth_week'];landmark=8
eligible=birth>=28;rows=np.arange(n)*17+landmark
dates=w[0]['dates'][c['entry_day']];year=dates.astype('datetime64[Y]').astype(int)+1970
# No future exposures or outcome-generated covariates enter landmark prediction.
X=np.column_stack([c['age_c'],c['smoking'],c['deprivation'],h['temp0'][rows],h['temp1'][rows],h['base'][rows,6:10]])
y=(birth<37).astype(int)
early=(dates+np.timedelta64(119,'D'))<=np.datetime64('2022-12-31')
late=dates>=np.datetime64('2023-01-01')
splits={'train':eligible&(c['site']<4)&early,'temporal':eligible&(c['site']<4)&late,'geographic':eligible&(c['site']>=4)&early,'combined_shift':eligible&(c['site']>=4)&late}
gap=eligible&~(early|late)
assert sum(s.sum() for s in splits.values())+gap.sum()==eligible.sum()
assert not np.any(splits['train']&splits['temporal'])
models={'Maternal baseline':make_pipeline(StandardScaler(),LogisticRegression(C=1,max_iter=2000,random_state=CONFIG['seed'])),'Logistic':make_pipeline(StandardScaler(),LogisticRegression(C=1,max_iter=2000,random_state=CONFIG['seed'])),'Random forest':RandomForestClassifier(n_estimators=400,min_samples_leaf=25,max_features=.7,n_jobs=1,random_state=CONFIG['seed']),'Gradient boosting':GradientBoostingClassifier(n_estimators=150,learning_rate=.03,max_depth=2,min_samples_leaf=25,random_state=CONFIG['seed'])}
metrics=[];predictions=[];rocs=[];calibration=[];rng=np.random.default_rng(CONFIG['seed']+77)
for name,model in models.items():
 features=X[:,:3] if name=='Maternal baseline' else X
 model.fit(features[splits['train']],y[splits['train']])
 for split in ['temporal','geographic','combined_shift']:
  sel=splits[split];yy=y[sel];pr=model.predict_proba(features[sel])[:,1];z=logit(np.clip(pr,1e-6,1-1e-6))
  cal=sm.GLM(yy,sm.add_constant(z),family=sm.families.Binomial()).fit();intercept=float(sm.GLM(yy,np.ones((len(yy),1)),offset=z,family=sm.families.Binomial()).fit().params[0])
  aucs=[]
  for b in range(300):
   j=rng.integers(0,len(yy),len(yy))
   if len(np.unique(yy[j]))==2:aucs.append(roc_auc_score(yy[j],pr[j]))
  lo,hi=np.quantile(aucs,[.025,.975])
  metrics.append(dict(model=name,test=split,n=len(yy),events=int(yy.sum()),prevalence=yy.mean(),auc=roc_auc_score(yy,pr),auc_lower=lo,auc_upper=hi,average_precision=average_precision_score(yy,pr),brier=brier_score_loss(yy,pr),calibration_intercept=intercept,calibration_slope=cal.params[1],mean_prediction=pr.mean()))
  for ident,yi,pi in zip(c['id'][sel],yy,pr):predictions.append(dict(id=ident,model=name,test=split,outcome=yi,prediction=pi))
  fpr,tpr,threshold=roc_curve(yy,pr)
  for a,b in zip(fpr,tpr):rocs.append(dict(model=name,test=split,fpr=a,tpr=b))
  bins=pd.qcut(pr,5,duplicates='drop')
  for _,g in pd.DataFrame(dict(pred=pr,y=yy,bin=bins)).groupby('bin',observed=True):
   calibration.append(dict(model=name,test=split,mean_prediction=g.pred.mean(),observed=g.y.mean(),n=len(g)))
pd.DataFrame(metrics).to_csv(ROOT/'results/prediction_metrics.csv',index=False);pd.DataFrame(predictions).to_csv(ROOT/'results/landmark_predictions.csv',index=False);pd.DataFrame(rocs).to_csv(ROOT/'results/roc.csv',index=False);pd.DataFrame(calibration).to_csv(ROOT/'results/calibration.csv',index=False)
pd.DataFrame([dict(split=k,n=int(v.sum()),events=int(y[v].sum())) for k,v in {**splits,'excluded_transition':gap}.items()]).to_csv(ROOT/'results/prediction_splits.csv',index=False)
# Independent R cross-check for an unpenalised logistic specification.
xx=np.column_stack([np.ones(n),X]);train=splits['train'];f=sm.GLM(y[train],xx[train],family=sm.families.Binomial()).fit(tol=1e-10)
pd.DataFrame(dict(term=np.arange(xx.shape[1]),coefficient=f.params)).to_csv(ROOT/'results/python_landmark_logistic.csv',index=False)
pd.DataFrame(np.column_stack([c['id'],train.astype(int),y,xx]),columns=['id','train','outcome']+[f'x{i}' for i in range(xx.shape[1])]).to_csv(ROOT/'data/landmark_design.csv',index=False)
# Fit the illustrative weekly model for explanatory response curves and hypothetical standardisation.
x,yy,_,keep=design(c,w,'nonlinear');fit=sm.GLM(yy,x,family=sm.families.Binomial(link=sm.families.links.CLogLog())).fit(tol=1e-10)
pd.DataFrame(dict(term=NAMES,coefficient=fit.params,se=fit.bse)).to_csv(ROOT/'results/illustrative_coefficients.csv',index=False)
curve=[]
for t in np.linspace(20,35,61):
 z=max(t-25,0)/5;contrast=np.r_[np.zeros(15),z,z*z,z,z*z];theta=contrast@fit.params;se=np.sqrt(contrast@fit.cov_params()@contrast)
 curve.append(dict(temperature=t,ratio=np.exp(theta),lower=np.exp(theta-1.96*se),upper=np.exp(theta+1.96*se),truth=np.exp(contrast[15:]@BETA)))
pd.DataFrame(curve).to_csv(ROOT/'results/response_curve.csv',index=False)
baseline=np.array([-5.62,.84,0,.12,.4,.25,-.025,.15,-.1,.025,.1,-.1,.05,-.05,.12])
risks=[];matrices=[]
for cooling in [0,2]:
 full=np.column_stack([h['base'],heat_basis(h['temp0']-cooling,h['temp1']-cooling)]);matrices.append(full)
 risk=1-np.prod(1-cloglog(full@fit.params).reshape(n,17),axis=1)
 truth=1-np.prod(1-cloglog(full@np.r_[baseline,BETA]).reshape(n,17),axis=1)
 risks.append((risk,truth))
gradients=[]
for matrix in matrices:
 e=np.exp(matrix@fit.params);H=e.reshape(n,17).sum(axis=1)
 gradients.append(((e[:,None]*matrix).reshape(n,17,matrix.shape[1]).sum(axis=1)*np.exp(-H)[:,None]).mean(axis=0))
gradient=gradients[1]-gradients[0];delta_se=float(np.sqrt(gradient@fit.cov_params()@gradient));delta=float((risks[1][0]-risks[0][0]).mean())
summary=dict(baseline_risk=float(risks[0][0].mean()),cooling_risk=float(risks[1][0].mean()),risk_difference=delta,truth_baseline_risk=float(risks[0][1].mean()),truth_cooling_risk=float(risks[1][1].mean()),truth_risk_difference=float((risks[1][1]-risks[0][1]).mean()),lower=delta-1.96*delta_se,upper=delta+1.96*delta_se,delta_method_se=delta_se,interval_method='Analytic delta method conditional on simulated population and weather',eligible=int(eligible.sum()),total=n,preterm_events=int((birth<37).sum()),scope='Hypothetical two-degree exposure reduction; full latent schedules used only for explanatory simulation standardisation, not prospective prediction')
(ROOT/'results/standardisation.json').write_text(json.dumps(summary,indent=2))
(ROOT/'provenance/prediction_specification.json').write_text(json.dumps(dict(features=['age_c','smoking','deprivation','T0','T1','rain','season_sin','season_cos','calendar_year'],training_only=True,tuning='Prespecified settings; no test-set tuning',bootstrap_auc=300,sklearn=metadata.version('scikit-learn'),model_parameters={k:str(m) for k,m in models.items()}),indent=2))
print(pd.DataFrame(metrics).to_string(index=False),flush=True)
