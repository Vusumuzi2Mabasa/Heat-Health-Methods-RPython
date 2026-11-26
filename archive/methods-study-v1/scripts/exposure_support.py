from common import *
import pandas as pd
w=weather();c=load_cohort(ROOT/'data/illustrative_cohort.csv');h=histories(c,w);p=h['person'];at_risk=h['interval']+20<=np.minimum(c['birth_week'][p],36)
rows=[]
for site in range(6):
 sel=at_risk&(c['site'][p]==site);t0=h['temp0'][sel];t1=h['temp1'][sel]
 rows.append(dict(site=w[site]['site_name'],intervals=int(sel.sum()),min_T0=t0.min(),max_T0=t0.max(),median_T0=np.median(t0),fraction_both_above25=float(((t0>25)&(t1>25)).mean()),fraction_both_atleast30=float(((t0>=30)&(t1>=30)).mean())))
pd.DataFrame(rows).to_csv(ROOT/'results/exposure_support.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
