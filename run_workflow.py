"""Compatibility entry point for the primary R association analysis."""
from pathlib import Path
import argparse,os,shutil,subprocess
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--rscript',default=shutil.which('Rscript'));p.add_argument('--quarto',default=shutil.which('quarto'));p.add_argument('--render',action='store_true');p.add_argument('--analysis-only',action='store_true');p.add_argument('--smoke',action='store_true',help='Compatibility alias for --analysis-only');a=p.parse_args()
if not a.rscript:p.error('Supply --rscript /path/to/Rscript or use Rscript --vanilla run_association.R')
env=os.environ.copy()
if a.quarto:env['QUARTO_PATH']=a.quarto
cmd=[a.rscript,'--vanilla',str(ROOT/'run_association.R')]
if a.analysis_only or a.smoke:cmd.append('--analysis-only')
if a.render:cmd.append('--render')
subprocess.run(cmd,cwd=ROOT,env=env,check=True)
