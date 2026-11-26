"""Audit exported scientific records and write reproducible content checksums."""
from common import *
import pandas as pd
import hashlib,platform
r=pd.read_csv(ROOT/'results/r_fits.csv',keep_default_na=False);py=pd.read_csv(ROOT/'results/python_fits.csv',keep_default_na=False)
keys=['replicate','scenario','model'];z=r.merge(py,on=keys,suffixes=('_R','_Python'),validate='one_to_one')
assert len(z)==1500
assert set(r.scenario)==set(py.scenario)==set(CONFIG['scenarios'])
assert (r.groupby(['scenario','model']).size()==100).all()
for field in ['events','rows','observed_events','truth']:
 np.testing.assert_allclose(z[field+'_R'],z[field+'_Python'],rtol=0,atol=1e-12)
assert r.error.eq('').all() and py.error.eq('').all()
for model in ['r_prediction_metrics.csv','prediction_metrics.csv']:
 d=pd.read_csv(ROOT/'results'/model);assert len(d)==12;assert d.auc.between(0,1).all();assert d.brier.between(0,1).all()
assert not any(x in (ROOT/'README.md').read_text() for x in ['ABSTRACT_RESULTS','SIMULATION_TABLE','FIGURES/','PREDICTION_NARRATIVE'])
release=dict(passed=True,matched_records=1500,scenario_labels_preserved=True,record_counts_identical=True,main_fits_without_errors=True,prediction_records_per_language=12,synthetic_health_data_only=True,python=platform.python_version())
(ROOT/'results/release_audit.json').write_text(json.dumps(release,indent=2))
files=[]
for folder in ['data','results','figures','scripts','provenance','report']:
 for path in sorted((ROOT/folder).glob('*')):
  if path.is_file() and not any(t in path.name for t in ['pilot_','smoke_','checksums','Manuscript.pdf']):files.append(path)
files+=[ROOT/'README.md',ROOT/'config.json',ROOT/'requirements.txt',ROOT/'run_workflow.py']
(ROOT/'provenance/checksums.sha256').write_text('\n'.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(ROOT).as_posix() for p in files)+'\n')
print(json.dumps(release,indent=2))
