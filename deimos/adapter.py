"""Release FITS -> project time/flux/error/quarter contract. No search or pixel extraction."""
import pathlib,json,re,time,hashlib
import numpy as np
import pandas as pd
from astropy.io import fits

def adapt(path, outdir=None, quarters=None, family='kbonus-bkg', quality_zero=True):
 start=time.perf_counter();path=pathlib.Path(path); frames=[];audit=[]
 with fits.open(path,memmap=False) as hdus:
  primary=dict(hdus[0].header);
  if primary.get('HLSPID','').strip().upper()!='KBONUS-BKG':raise ValueError('wrong archive family')
  headers=[x.header.tostring(sep='\n',endcard=True,padding=False) for x in hdus]
  for i,h in enumerate(hdus[1:],1):
   match=re.fullmatch(r'LIGHTCURVE_Q(\d+)',h.name)
   if family=='kbonus-bkg' and not match:continue
   if family=='kbonus-apexba' and not hasattr(h,'columns'):continue
   q=int(match[1]) if match else int(h.header.get('QUARTER',primary.get('QUARTER',-1)))
   if quarters is not None and q not in quarters:continue
   columns={c.upper():c for c in h.columns.names};needed=['TIME','FLUX','FLUX_ERR']
   if not set(needed)<=columns.keys():raise ValueError('missing adapter columns')
   if any(str(h.columns[columns[c]].unit).strip()!='e-/s' for c in ['FLUX','FLUX_ERR']):raise ValueError('unverified flux units')
   hdr=dict(h.header);ref=float(hdr.get('BJDREFI',primary.get('BJDREFI',2454833)))+float(hdr.get('BJDREFF',primary.get('BJDREFF',0)))
   timezero=float(hdr.get('TIMEZERO',primary.get('TIMEZERO',0)));ts=hdr.get('TIMESYS',primary.get('TIMESYS'))
   if not ts or ts.strip().upper()!='TDB':raise ValueError('not documented TDB '+str(ts))
   unit=str(h.columns[columns['TIME']].unit).strip()
   if unit!='BJD - 2454833':raise ValueError('unverified TIME unit '+unit)
   # Archive TIME is relative days. Accept only verified 2454833 reference; do not guess other offsets.
   if ref!=2454833:raise ValueError('unexpected BJD reference '+str(ref))
   raw={c:np.array(h.data[c]) for c in h.columns.names if np.array(h.data[c]).ndim==1}
   frame=pd.DataFrame(raw).add_prefix('original_')
   t=np.array(h.data[columns['TIME']],float)+timezero
   y=np.array(h.data[columns['FLUX']],float);dy=np.array(h.data[columns['FLUX_ERR']],float)
   qualitycol=columns.get('SAP_QUALITY') if family=='kbonus-bkg' else columns.get('QUALITY')
   quality=np.array(h.data[qualitycol],dtype=np.int64) if qualitycol else np.full(len(t),-1)
   valid=np.isfinite(t)&np.isfinite(y)&np.isfinite(dy)&(dy>0)
   if quality_zero:valid&=quality==0
   frac=float(hdr.get('PSFFRAC',np.nan))
   if family=='kbonus-bkg' and (not np.isfinite(frac) or frac<=.5):valid[:]=False
   norm=float(np.median(y[valid])) if valid.any() else np.nan
   if not np.isfinite(norm) or norm<=0:valid[:]=False
   frame['time']=t;frame['flux']=y/norm if np.isfinite(norm) and norm>0 else np.nan;frame['flux_err']=dy/norm if np.isfinite(norm) and norm>0 else np.nan
   frame['quarter']=q;frame['quality']=quality;frame['adapter_valid']=valid;frame['normalization_median_e_s']=norm
   frame['input_sha256']=hashlib.sha256(path.read_bytes()).hexdigest();frame['input_hdu']=i;frame['input_row']=np.arange(len(frame));frames.append(frame)
   audit.append({'hdu':i,'quarter':q,'rows':len(frame),'valid_rows':int(valid.sum()),'time_min_bkjd':float(np.nanmin(t)) if len(t) else None,'time_max_bkjd':float(np.nanmax(t)) if len(t) else None,'median_cadence_days':float(np.nanmedian(np.diff(np.sort(t[np.isfinite(t)])))),'median_flux_e_s':norm,'quality_column':qualitycol,'quality_counts':{str(k):int(v) for k,v in zip(*np.unique(quality,return_counts=True))},'finite_flux_rows':int(np.sum(np.isfinite(y))),'nonpositive_flux_rows':int(np.sum(np.isfinite(y)&(y<=0))),'finite_positive_error_rows':int(np.sum(np.isfinite(dy)&(dy>0))),'header':{k:hdr[k] for k in hdr if k not in ('COMMENT','HISTORY')},'columns':[{'name':c.name,'format':c.format,'unit':c.unit} for c in h.columns]})
 if not frames:raise ValueError('No supported quarter extensions')
 frame=pd.concat(frames,ignore_index=True)
 kept=frame[frame.adapter_valid].sort_values(['time','quarter','input_row']).drop_duplicates('time').copy()
 if len(kept)==0:raise ValueError('No usable cadences')
 meta={'family':family,'input':str(path),'input_bytes':path.stat().st_size,'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'original_header_cards':headers,'primary_header':{k:v for k,v in primary.items() if k not in ('COMMENT','HISTORY')},'quarter_audit':audit,'preserved_rows':len(frame),'normalized_valid_rows':len(kept),'adapter_seconds':time.perf_counter()-start,'time_system':'BKJD = BJD_TDB - 2454833 days','flux_units':'relative quarter median','error_units':'relative; formal archive uncertainties, covariance not supplied','quality_policy':'strict zero original TPF flags; deliberately conservative engineering policy, not PDCSAP default','upstream_dependency':'corrected PSF model fit to quarter-wide pixels and shared TPF stack; not independent of upstream processing','flattened_stitched_used':False}
 if outdir:
  outdir=pathlib.Path(outdir);outdir.mkdir(parents=True,exist_ok=True);stem=path.name.removesuffix('.fits')
  frame.to_csv(outdir/(stem+'-preserved.csv.gz'),index=False);kept.to_csv(outdir/(stem+'-normalized.csv.gz'),index=False);(outdir/(stem+'-provenance.json')).write_text(json.dumps(meta,indent=2,default=str,allow_nan=True));meta['adapter_and_export_seconds']=time.perf_counter()-start
 return kept,meta

if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('--out',required=True);p.add_argument('--family',default='kbonus-bkg');a=p.parse_args();_,m=adapt(a.input,a.out,family=a.family);print(json.dumps({k:m[k] for k in ['input_bytes','normalized_valid_rows','adapter_seconds']},indent=2))
