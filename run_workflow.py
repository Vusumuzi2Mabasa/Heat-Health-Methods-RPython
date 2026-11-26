"""Rebuild a simulation smoke test or the complete manuscript study."""
from pathlib import Path
import os,sys,subprocess,argparse,shutil,json
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--smoke',action='store_true');p.add_argument('--download-weather',action='store_true');p.add_argument('--rscript',default=shutil.which('Rscript'));p.add_argument('--quarto',default=shutil.which('quarto'));p.add_argument('--render',action='store_true');a=p.parse_args()
if not a.rscript:sys.exit('Supply --rscript /path/to/Rscript or add Rscript to PATH')
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',R_ENVIRON_USER='NUL' if os.name=='nt' else '/dev/null',R_PROFILE_USER='NUL' if os.name=='nt' else '/dev/null')
def run(*cmd):print('Running:',*cmd,flush=True);subprocess.run([str(x) for x in cmd],cwd=ROOT,check=True)
def py(name,*args):run(sys.executable,ROOT/'scripts'/name,*args)
if a.download_weather:py('fetch_weather.py')
if a.smoke:
 py('generate.py','--reps','2','--n','1800','--out','cache/smoke','--manifest','data/smoke_manifest.csv')
 py('analyse.py','--manifest','data/smoke_manifest.csv','--output','results/smoke_python.csv')
 run(a.rscript,'--vanilla',ROOT/'scripts/analyse.R',ROOT,'data/smoke_manifest.csv','results/smoke_r.csv')
 py('check_smoke.py');sys.exit(0)
py('generate.py');py('analyse.py');run(a.rscript,'--vanilla',ROOT/'scripts/analyse.R',ROOT)
py('prediction.py');run(a.rscript,'--vanilla',ROOT/'scripts/prediction_check.R',ROOT);run(a.rscript,'--vanilla',ROOT/'scripts/prediction_models.R',ROOT)
py('exposure_support.py');py('test_invariants.py');py('summarise.py');run(a.rscript,'--vanilla',ROOT/'scripts/figures.R',ROOT);run(a.rscript,'--vanilla',ROOT/'scripts/environment.R',ROOT);py('write_manuscript.py')
if a.render:
 if not a.quarto:sys.exit('Supply --quarto /path/to/quarto for Word rendering')
 run(a.quarto,'render',ROOT/'report/Heat_Health_Manuscript.qmd','--to','docx')
 py('style_report.py',ROOT/'report/Heat_Health_Manuscript.docx')
 py('release_audit.py')
print('Finished. DOCX visual QA is a separate release step.')
