from common import *
import argparse,time,platform
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm
parser=argparse.ArgumentParser();parser.add_argument('--manifest',default='data/simulation_manifest.csv');parser.add_argument('--output',default='results/python_fits.csv');a=parser.parse_args()
w=weather();manifest=pd.read_csv(ROOT/a.manifest,keep_default_na=False);records=[];coeff=[];started=time.time()
for _,row in manifest.iterrows():
 c=load_cohort(ROOT/row.path)
 cached=None
 for model in ['nonlinear','linear','oracle']:
  x,y,h,keep=design(c,w,model)
  try:
   fit=cached if model=='oracle' and row.scenario in ['null','moderate','nonlinear'] else sm.GLM(y,x,family=sm.families.Binomial(link=sm.families.links.CLogLog())).fit(maxiter=100,tol=1e-10)
   if model=='nonlinear':cached=fit
   contrast=np.r_[np.zeros(15),np.ones(x.shape[1]-15)];theta=float(contrast@fit.params);se=float(np.sqrt(contrast@fit.cov_params()@contrast))
   records.append(dict(replicate=row.replicate,scenario=row.scenario,model=model,n=row.n,events=row.events,rows=len(y),observed_events=int(y.sum()),truth=row.truth,estimate=theta,se=se,lower=theta-1.96*se,upper=theta+1.96*se,p=2*norm.sf(abs(theta/se)),converged=bool(fit.converged),loglik=fit.llf,error=''))
   if row.replicate==1:
    for j,b in enumerate(fit.params):coeff.append(dict(scenario=row.scenario,model=model,term=j,estimate=b,se=fit.bse[j]))
  except Exception as e:
   records.append(dict(replicate=row.replicate,scenario=row.scenario,model=model,converged=False,error=str(e)))
 if int(row.replicate)%10==0 and row.scenario==manifest.scenario.iloc[-1]:print('Python replicate',row.replicate,'seconds',round(time.time()-started),flush=True)
out=ROOT/a.output;out.parent.mkdir(parents=True,exist_ok=True);pd.DataFrame(records).to_csv(out,index=False);pd.DataFrame(coeff).to_csv(out.with_name(out.stem+'_coefficients.csv'),index=False)
(ROOT/'provenance'/f'{out.stem}_runtime.json').write_text(json.dumps(dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,statsmodels=sm.__version__,manifest=a.manifest,output=a.output,fit_records=len(records),elapsed_seconds=time.time()-started),indent=2))
print('Python fits',len(records),'seconds',round(time.time()-started),flush=True)
