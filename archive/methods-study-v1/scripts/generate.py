from common import *
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--reps',type=int,default=CONFIG['replicates']);parser.add_argument('--n',type=int,default=CONFIG['simulation_n']);parser.add_argument('--out',default='cache/cohorts');parser.add_argument('--manifest',default='data/simulation_manifest.csv');a=parser.parse_args()
out=ROOT/a.out;out.mkdir(parents=True,exist_ok=True);w=weather()
low=int((np.datetime64('2020-01-15')-w[0]['dates'][0]).astype(int));high=int((np.datetime64('2024-07-01')-w[0]['dates'][0]).astype(int))
def cohort(seed,n,scenario):
 rng=np.random.default_rng(seed);site=np.tile(np.arange(6),int(np.ceil(n/6)))[:n];rng.shuffle(site)
 c=dict(id=np.arange(n),site=site,entry_day=rng.integers(low,high,n),age_c=(np.clip(rng.normal(29,6,n),18,44)-30)/10,smoking=rng.binomial(1,.15,n),deprivation=rng.binomial(1,.4,n),personal_offset=np.zeros(n),birth_week=np.full(n,37),missing_mask=np.zeros(n,dtype=np.int64))
 offsets=rng.normal(0,2,n);uniform=rng.random((n,17));missing_u=rng.random((n,17))
 if scenario=='exposure_error':c['personal_offset']=offsets
 h=histories(c,w);base=h['base'];p=h['person'];b=BETA*(0 if scenario=='null' else .5 if scenario=='moderate' else 1)
 target=np.array([-5.62,.84,0,.12,.4,.25,-.025,.15,-.1,.025,.1,-.1,.05,-.05,.12])
 heat=heat_basis(h['temp0']+c['personal_offset'][p],h['temp1']+c['personal_offset'][p]);eta=base@target+heat@b;prob=cloglog(eta).reshape(n,17)
 occurred=uniform<prob;first=np.argmax(occurred,axis=1);has=occurred.any(axis=1);c['birth_week'][has]=20+first[has]
 if scenario=='informative_missingness':
  event=(np.arange(20,37)[None,:]==c['birth_week'][:,None]);z0=heat[:,0].reshape(n,17)
  pm=1/(1+np.exp(-(-3.2+.45*z0+.5*c['deprivation'][:,None]+1.5*event)))
  missing=missing_u<pm;c['missing_mask']=(missing*(1<<np.arange(17))).sum(axis=1)
 return c,float(b.sum())
def write(path,c):
 opener=gzip.open if str(path).endswith('.gz') else open
 with opener(path,'wt',newline='') as f:
  wr=csv.writer(f);wr.writerow(c.keys());wr.writerows(zip(*c.values()))
manifest=[]
for rep in range(1,a.reps+1):
 for scenario in CONFIG['scenarios']:
  c,truth=cohort(CONFIG['seed']+rep,a.n,scenario);p=out/f'{rep:03d}_{scenario}.csv.gz';write(p,c)
  manifest.append(dict(replicate=rep,scenario=scenario,path=p.relative_to(ROOT).as_posix(),n=a.n,events=int((c['birth_week']<37).sum()),truth=truth,seed=CONFIG['seed']+rep))
 if rep%10==0:print('Generated replicate',rep,flush=True)
mp=ROOT/a.manifest;mp.parent.mkdir(parents=True,exist_ok=True)
with mp.open('w',newline='') as f:wr=csv.DictWriter(f,fieldnames=manifest[0]);wr.writeheader();wr.writerows(manifest)
if a.reps>=CONFIG['replicates']:
 c,truth=cohort(CONFIG['seed']+9000,CONFIG['illustrative_n'],'nonlinear');write(ROOT/'data/illustrative_cohort.csv',c)
json_path=ROOT/'provenance/generating_parameters.json';json_path.write_text(json.dumps(dict(config=CONFIG,baseline_coefficients=[-5.62,.84,0,.12,.4,.25,-.025,.15,-.1,.025,.1,-.1,.05,-.05,.12],baseline_names=NAMES[:15],weekly_hazard_link='complementary log-log',weather_real=True,health_data_simulated=True),indent=2))
print('Saved cohorts:',len(manifest),flush=True)
