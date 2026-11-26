from common import *
import pandas as pd
a=pd.read_csv(ROOT/'results/smoke_python.csv',keep_default_na=False);b=pd.read_csv(ROOT/'results/smoke_r.csv',keep_default_na=False)
assert len(a)==len(b)==30 and a.converged.all() and b.converged.all()
assert set(a.scenario)==set(b.scenario)==set(CONFIG['scenarios'])
checks={}
for k,tol in [('estimate',1e-5),('se',1e-4),('loglik',1e-5)]:
 diff=float((a[k]-b[k]).abs().max());assert diff<tol,(k,diff);checks[k]=diff;print(k,'maximum absolute difference',diff)
(ROOT/'results/smoke_validation.json').write_text(json.dumps(dict(passed=True,records_per_language=30,scenarios=CONFIG['scenarios'],maximum_absolute_differences=checks),indent=2))
print('SMOKE TEST PASSED: 30 matched fit records per language')
