from common import *
import unittest
class Invariants(unittest.TestCase):
 def test_scenario_names(self):
  import pandas as pd
  d=pd.read_csv(ROOT/'data/simulation_manifest.csv',keep_default_na=False)
  self.assertEqual(set(d.scenario),set(CONFIG['scenarios']));self.assertEqual(int((d.scenario=='null').sum()),100)
 def test_weather(self):
  w=weather();self.assertEqual(len(w),6)
  for z in w:
   self.assertTrue(np.isfinite(z['tx']).all());self.assertTrue((z['rain']>=0).all())
 def test_windows(self):
  c=load_cohort(ROOT/'data/illustrative_cohort.csv');w=weather();h=histories(c,w)
  for i in [0,170,9000]:
   p=h['person'][i];d=h['day'][i];z=w[c['site'][p]]
   self.assertAlmostEqual(h['temp0'][i],z['tx'][d-7:d].mean(),places=9)
   self.assertAlmostEqual(h['temp1'][i],z['tx'][d-14:d-7].mean(),places=9)
 def test_risk_sets(self):
  c=load_cohort(ROOT/'data/illustrative_cohort.csv');x,y,h,k=design(c,weather());self.assertEqual(int(y.sum()),int((c['birth_week']<37).sum()))
  self.assertTrue(np.all(h['interval'][k]+20<=c['birth_week'][h['person'][k]]))
 def test_contrast(self):
  np.testing.assert_allclose(heat_basis(np.array([25,30]),np.array([25,30])),[[0,0,0,0],[1,1,1,1]])
  self.assertAlmostEqual(float(BETA.sum()),.42)
 def test_landmark(self):
  import pandas as pd
  d=pd.read_csv(ROOT/'data/landmark_design.csv');c=load_cohort(ROOT/'data/illustrative_cohort.csv');s=pd.read_csv(ROOT/'results/prediction_splits.csv')
  self.assertEqual(int(s.n.sum()),int((c['birth_week']>=28).sum()))
  train=d.train.to_numpy()==1;dates=np.datetime64('2019-01-01')+c['entry_day'].astype('timedelta64[D]')
  self.assertTrue(np.all(dates[train]+np.timedelta64(119,'D')<=np.datetime64('2022-12-31')))
  self.assertTrue(np.all(c['birth_week'][train]>=28));self.assertTrue(np.all(c['site'][train]<4))
 def test_missingness_oracle(self):
  c=load_cohort(ROOT/'cache/cohorts/001_informative_missingness.csv.gz');w=weather();xo,yo,_,_=design(c,w,'oracle');xm,ym,_,_=design(c,w,'nonlinear')
  self.assertEqual(int(yo.sum()),int((c['birth_week']<37).sum()));self.assertLess(len(ym),len(yo));self.assertLess(ym.sum(),yo.sum())
if __name__=='__main__':unittest.main()
