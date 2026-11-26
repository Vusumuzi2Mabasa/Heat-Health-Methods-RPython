from pathlib import Path
from urllib.request import urlopen, Request
from urllib.parse import urlencode
import csv,json,hashlib,time
root=Path(__file__).resolve().parents[1];config=json.loads((root/'config.json').read_text())
sites=[('Johannesburg',-26.2041,28.0473,'Africa/Johannesburg'),('Durban',-29.8587,31.0218,'Africa/Johannesburg'),('Cape_Town',-33.9249,18.4241,'Africa/Johannesburg'),('Harare',-17.8252,31.0335,'Africa/Harare'),('Nairobi',-1.2921,36.8219,'Africa/Nairobi'),('Abidjan',5.36,-4.0083,'Africa/Abidjan')]
cache=root/'cache/weather';cache.mkdir(parents=True,exist_ok=True)
data=root/'data';data.mkdir(exist_ok=True);prov=root/'provenance';prov.mkdir(exist_ok=True)
records=[];manifest=[]
for i,(name,lat,lon,tz) in enumerate(sites):
 params=dict(latitude=lat,longitude=lon,start_date=config['weather_start'],end_date=config['weather_end'],daily='temperature_2m_max,temperature_2m_min,precipitation_sum',models=config['climate_model'],timezone=tz)
 url='https://archive-api.open-meteo.com/v1/archive?'+urlencode(params);p=cache/(name+'_'+config['climate_model']+'.json')
 if not p.exists():
  for attempt in range(4):
   try:
    with urlopen(Request(url,headers={'User-Agent':'Heat-Health-Methods-RPython research portfolio'}),timeout=120) as response:raw=response.read()
    payload=json.loads(raw)
    if 'daily' not in payload:raise RuntimeError(str(payload))
    p.write_bytes(raw);break
   except Exception:
    if attempt==3:raise
    time.sleep(2*(attempt+1))
 payload=json.loads(p.read_bytes());daily=payload['daily']
 assert len(daily['time'])==2192
 for k,date in enumerate(daily['time']):
  mx,mn,rain=[daily[x][k] for x in ['temperature_2m_max','temperature_2m_min','precipitation_sum']]
  assert None not in [mx,mn,rain] and -30<=mn<=mx<=55 and rain>=0
  records.append(dict(site=i,site_name=name,date=date,tmax_c=mx,tmin_c=mn,rain_mm=rain,latitude=lat,longitude=lon))
 manifest.append(dict(site=name,url=url,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size,latitude_returned=payload['latitude'],longitude_returned=payload['longitude'],timezone=payload['timezone']))
 print('Verified weather:',name,flush=True)
with (data/'weather.csv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=records[0]);writer.writeheader();writer.writerows(records)
(prov/'weather_manifest.json').write_text(json.dumps({'provider':'Open-Meteo Historical Weather API','model_requested':config['climate_model'],'weather_type':'Gridded reanalysis estimates, not personal exposure or station observations','license':'CC BY 4.0; retain Open-Meteo and underlying ERA5 attribution','source_documentation':'https://open-meteo.com/en/docs/historical-weather-api','sites_are':'Illustrative city coordinates, not Wits PHR recruitment facilities','files':manifest},indent=2))
print('Weather rows:',len(records),flush=True)
