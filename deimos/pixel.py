"""Pixel host test: which star in the Kepler pixels does a transit dip come from?

Usage: python -m deimos.pixel LEADS_CSV --out DIR [RANK ...]
LEADS_CSV columns: rank,lg,role,gaia,origin,peak_index,src (src = folder with that source's
real.json, receipt.json, split.json). Downloads are cached under DIR.
"""
import argparse, hashlib, json, os, pathlib, traceback, urllib.parse, urllib.request, warnings
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from scipy.interpolate import RectBivariateSpline
from deimos.core import digest, fixed_holdout, serial, transit_masks

# Neighbours whose PRF lies mostly off the stamp have unrecoverable injected dips, so the
# extraction control would fail. Before localize_tpf, a neighbour is kept as a hypothesis only if
# >= 0.25 of its PRF falls on the stamp (outcome-blind: geometry only). Dropped stars are listed
# with G and in-stamp fraction.
#
# Decision rule (fixed before any result is seen):
#  A. Difference image, per quarter: a star is assigned if localize_tpf reports
#     *_supported or beats the next hypothesis by delta chi2 >= 9.
#     Lead-level: assigned star = the one assigned in >= 2/3 of assigned quarters,
#     with no quarter assigning a different star.
#  B. Extraction, only in quarters whose source-injection control passed: per-source
#     fractional dip at the frozen ephemeris, z against 24 off-transit phase offsets.
#     Assigned star = the only source with z >= 5 while every other source has z < 3.
#  KEEP if A and B both give the target. DROP if both give the same other star.
#  Otherwise UNRESOLVED. No refit, no period search, no holdout quarters.

def source_hash(path):
    digest_ = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest_.update(block)
    return digest_.hexdigest()


def json_write(path, value) -> None:
    """Atomically replace a JSON artifact."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, default=str, allow_nan=False))
    tmp.replace(path)


def csv_write(path, frame) -> None:
    """Atomically replace a CSV artifact."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.csv.tmp')
    frame.to_csv(tmp, index=False)
    tmp.replace(path)

def read(path,default=None):
    return json.loads(path.read_text()) if path.exists() else default


def freeze(path,value):
    value=serial(value)
    if path.exists():
        if digest(read(path))!=digest(value): raise ValueError(f'Frozen artifact drift: {path.name}')
    else: json_write(path,value)
    return value

def official_prf(path, column, row, module=None, output=None):
    with fits.open(path) as hdus:
        if module is not None and hdus[0].header['MODULE'] != module: raise ValueError('PRF module mismatch')
        if output is not None and hdus[0].header['OUTPUT'] != output: raise ValueError('PRF output mismatch')
        images=[]; weights=[]
        for h in hdus[1:6]:
            images.append(np.asarray(h.data, float))
            weights.append(1/max(1e-6,np.hypot(column-h.header['CRVAL1P'],row-h.header['CRVAL2P'])))
        dx,dy = hdus[1].header['CDELT1P'],hdus[1].header['CDELT2P']
    image=np.average(images,axis=0,weights=weights)
    image/=image.sum()*dx*dy
    x=(np.arange(.5,image.shape[1]+.5)-image.shape[1]/2)*dx
    y=(np.arange(.5,image.shape[0]+.5)-image.shape[0]/2)*dy
    spline=RectBivariateSpline(y,x,image)
    def template(shape, origin, position):
        # Absolute pixel edge origin; WCS positions are converted to centers +0.5.
        col0,row0=origin; col,rr=position
        # The calibrated PRF already predicts a pixel response. Integrating it
        # over each pixel again broadens the source (it is not an optical PSF).
        xx = col0 + np.arange(shape[1]) + .5 - col
        yy = row0 + np.arange(shape[0]) + .5 - rr
        value = spline(yy, xx)
        value[(np.abs(yy) > max(abs(y)))[:, None] | (np.abs(xx) > max(abs(x)))[None, :]] = 0
        return value
    return template


def propagate_gaia(record, bkjd):
    pm=(record.get('pmra'),record.get('pmdec'),record.get('ref_epoch'))
    if any(v is None or not np.isfinite(float(v)) for v in pm):
        return float(record['ra']),float(record['dec']), 'missing_proper_motion_unresolved'
    coord=SkyCoord(ra=float(record['ra'])*u.deg,dec=float(record['dec'])*u.deg,
                   pm_ra_cosdec=float(pm[0])*u.mas/u.yr,pm_dec=float(pm[1])*u.mas/u.yr,
                   obstime=Time(float(pm[2]),format='jyear'))
    moved=coord.apply_space_motion(new_obstime=Time(float(bkjd)+2454833,format='jd',scale='tdb'))
    return moved.ra.deg,moved.dec.deg,'proper_motion_propagated_no_radial_velocity'


def compare_hypotheses(image,error,templates,astrometry_complete=False):
    valid=np.isfinite(image)&np.isfinite(error)&(error>0)
    if valid.sum()<8 or len(templates)<2: return {'status':'unresolved','reason':'Insufficient pixels or alternative hypotheses'}
    y=image[valid]; e=error[valid]; rows=[]
    for name,template in templates.items():
        x=np.column_stack([template[valid],np.ones(valid.sum())])/e[:,None]
        if np.linalg.matrix_rank(x)<2: continue
        cov=np.linalg.inv(x.T@x); coef=cov@x.T@(y/e)
        residual=(y-(np.column_stack([template[valid],np.ones(valid.sum())])@coef))/e
        rows.append({'hypothesis':name,'amplitude':coef[0],'amplitude_error':np.sqrt(cov[0,0]),
                     'snr':coef[0]/np.sqrt(cov[0,0]),'chi2':sum(residual**2),
                     'reduced_chi2':sum(residual**2)/(valid.sum()-2),'residual_rms_sigma':np.sqrt(np.mean(residual**2))})
    rows=sorted(rows,key=lambda r:r['chi2'])
    if len(rows)<2: return {'status':'unresolved','reason':'Degenerate source templates'}
    delta=rows[1]['chi2']-rows[0]['chi2']; best=rows[0]
    strong=bool(astrometry_complete and delta>=25 and best['snr']>=7 and best['reduced_chi2']<=3)
    status='target_supported' if strong and best['hypothesis']=='target' else 'off_target_supported' if strong else 'unresolved'
    return serial({'status':status,'fits':rows,'delta_chi2_next':delta,'astrometry_complete':astrometry_complete,
                   'valid_pixels':int(valid.sum()),'crowded':len(templates)>1,
                   'interpretation':'Fixed source hypotheses; PRF mismatch and correlated pixel errors limit significance. No planet confirmation.'})


def extract_sources(cube,errors,templates):
    """Optional linear PRF photometry; missing/degenerate sources fail closed."""
    design=np.column_stack([t.ravel() for t in templates.values()]+[np.ones(cube.shape[1]*cube.shape[2])])
    estimates=[]; uncertainties=[]
    for image,error in zip(cube,errors):
        good=np.isfinite(image.ravel())&np.isfinite(error.ravel())&(error.ravel()>0)
        if good.sum()<=design.shape[1]: raise ValueError('Insufficient source-extraction pixels')
        x=design[good]/error.ravel()[good,None]; y=image.ravel()[good]/error.ravel()[good]
        if np.linalg.cond(x)>1e6: raise ValueError('Degenerate source confusion')
        cov=np.linalg.inv(x.T@x); beta=cov@x.T@y
        estimates.append(beta[:-1]); uncertainties.append(np.sqrt(np.diag(cov))[:-1])
    return np.array(estimates),np.array(uncertainties)


def extraction_control(templates,seed=20261013):
    rng=np.random.default_rng(seed); shapes=list(templates.values()); n=160
    errors=np.full((n,*shapes[0].shape),.5)
    cube=np.zeros_like(errors)+10
    truth=np.ones((n,len(shapes)))*1000
    # Independent artificial dip assigned to each source in a separate time interval.
    for j in range(len(shapes)): truth[10+j*20:20+j*20,j]*=.98
    for j,t in enumerate(shapes): cube+=truth[:,j,None,None]*t
    cube+=rng.normal(0,errors)
    flux,err=extract_sources(cube,errors,templates)
    bias=np.median(flux-truth,axis=0)
    recovery=[]
    for j in range(len(shapes)):
        mask=np.zeros(n,bool);mask[10+j*20:20+j*20]=True
        recovered=np.median(flux[~mask,j])-np.median(flux[mask,j])
        recovery.append(float(recovered))
    passed=bool(np.all(np.abs(bias)<3*np.median(err,axis=0)) and all(abs(x-20)<10 for x in recovery))
    return serial({'passed':passed,'source_depths_recovered':recovery,'true_depth_flux':20,'bias':bias,
        'scope':'Synthetic source injection using official PRF templates; does not validate real background-star detrending or completeness.'})


def select_background_sources(neighbours,already_inspected=(),seed=2026101303,maximum=5):
    """Outcome-blind optional source selection; no dip/depth/period fields used."""
    eligible=[r for r in neighbours if str(r.get('source_id')) not in set(map(str,already_inspected))
        and r.get('phot_g_mean_mag') is not None and 11<=r['phot_g_mean_mag']<=16
        and 2<float(r.get('separation_arcsec',0))<=30 and r.get('pmra') is not None and r.get('pmdec') is not None]
    eligible=sorted(eligible,key=lambda r:str(r['source_id']))
    rng=np.random.default_rng(seed)
    return [eligible[i] for i in rng.permutation(len(eligible))[:maximum]]


def extract_background_discovery(cube,errors,times,templates,source_names,validation):
    """Validated optional extraction API; full real source validation required first."""
    if not validation.get('source_injection_passed') or not validation.get('known_source_timing_passed'):
        raise ValueError('Source extraction not validated for discovery use')
    if any(name not in templates for name in source_names): raise ValueError('Selected source absent from hypotheses')
    flux,error=extract_sources(cube,errors,templates)
    return {name:{'time':np.asarray(times),'flux':flux[:,list(templates).index(name)],
                   'flux_err':error[:,list(templates).index(name)],'source_specific':True} for name in source_names}

def difference_image(tpf_flux, time, candidate: dict, aperture_mask):
    """Build a locally matched per-event out-minus-in TPF image with per-pixel errors."""
    pix = np.asarray(tpf_flux, float)
    t = np.asarray(time, float)
    if pix.ndim != 3 or len(t) != len(pix):
        return {'status': 'manual_check_required', 'reason': 'TPF cube/time shape mismatch'}
    aperture = np.asarray(aperture_mask, bool)
    if aperture.shape != pix.shape[1:] or not aperture.any():
        aperture = np.ones(pix.shape[1:], dtype=bool)
    phase, events, inside = transit_masks(t, candidate['period_days'], candidate['epoch_bkjd'], candidate['duration_days'])
    event_differences, event_variances = [], []
    for event in np.unique(events[inside]):
        local_in = inside & (events == event)
        local_out = (events == event) & (np.abs(phase) >= candidate['duration_days'] / 2) & (np.abs(phase) < 2 * candidate['duration_days'])
        if local_in.sum() < 1 or local_out.sum() < 2:
            continue
        ins, outs = pix[local_in], pix[local_out]
        ni, no = np.isfinite(ins).sum(axis=0), np.isfinite(outs).sum(axis=0)
        with np.errstate(invalid='ignore', divide='ignore'), warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            mi = np.nansum(ins, axis=0) / ni
            mo = np.nansum(outs, axis=0) / no
            out_var = np.nanvar(outs, axis=0, ddof=1)
            in_var = np.nanvar(ins, axis=0, ddof=1)
            in_var = np.where(np.isfinite(in_var), in_var, out_var)
            vi = in_var / ni
            vo = out_var / no
        di = mo - mi
        var = vi + vo
        usable = (ni >= 1) & (no >= 2) & np.isfinite(di) & np.isfinite(var)
        di[~usable] = np.nan; var[~usable] = np.nan
        event_differences.append(di); event_variances.append(var)
    if len(event_differences) < 3:
        return {'status': 'manual_check_required', 'reason': 'Fewer than three events have local in/out baselines'}
    ds, vs = np.asarray(event_differences), np.asarray(event_variances)
    counts = np.isfinite(ds).sum(axis=0)
    with np.errstate(invalid='ignore', divide='ignore'), warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        diff = np.nansum(ds, axis=0) / counts
        within = np.nansum(vs, axis=0) / counts**2
        between = np.nanvar(ds, axis=0, ddof=1) / counts
    error = np.sqrt(within + np.nan_to_num(between, nan=0.0))
    invalid = counts < 2
    diff[invalid] = np.nan; error[invalid] = np.nan
    significance = diff / np.where(error > 0, error, np.nan)
    ap_valid = aperture & np.isfinite(significance)
    positive = ap_valid & (diff > 0) & (significance >= 3)
    return {'status': 'screened' if np.isfinite(significance).any() else 'manual_check_required',
            'difference_flux': diff, 'uncertainty': error, 'significance': significance,
            'events_used': len(event_differences), 'pixels_with_finite_estimate': int(np.isfinite(significance).sum()),
            'positive_significant_aperture_pixels': int(positive.sum()),
            'uncertainty_method': 'per-event local out/in mean errors plus event-to-event scatter; per pixel',
            'localization_interpretation': 'Qualitative pixel-level screen only; significant pixels do not prove the variable source is the target.'}

def localize_tpf(path,candidate,neighbours,out,limits):
    import lightkurve as lk
    tpf=lk.KeplerTargetPixelFile(path,quality_bitmask='default')
    t=np.asarray(tpf.time.bkjd);cube=np.asarray(tpf.flux.value,float);err=np.asarray(tpf.flux_err.value,float)
    hd=tpf.hdu[0].header;module,output=hd['MODULE'],hd['OUTPUT'];quarter=int(tpf.quarter)
    prfname=f'kplr{module:02d}.{output}_2011265_prf.fits'
    prfpath=limits.download('https://archive.stsci.edu/missions/kepler/fpc/prf/'+prfname,out/'prf'/prfname,30_000_000)
    prf=official_prf(prfpath,float(tpf.column)+cube.shape[2]/2,float(tpf.row)+cube.shape[1]/2,module,output)
    origin=(float(tpf.column),float(tpf.row));templates={};positions={}
    ra,dec=float(hd['RA_OBJ']),float(hd['DEC_OBJ'])
    host=next((n for n in neighbours if n.get('separation_arcsec',999)<2),None)
    astro=host is not None
    uncertainties={}
    if host:
        ra,dec,status=propagate_gaia(host,float(np.median(t)))
        if status.startswith('missing') or any(host.get(v) is None for v in ('ra_error','dec_error','pmra_error','pmdec_error')): astro=False
    x,y=tpf.wcs.all_world2pix([[ra,dec]],0)[0]
    positions['target']=(origin[0]+x+.5,origin[1]+y+.5)
    for i,n in enumerate(neighbours):
        if n.get('separation_arcsec',0)<2: continue
        nra,ndec,status=propagate_gaia(n,float(np.median(t)))
        x,y=tpf.wcs.all_world2pix([[nra,ndec]],0)[0]
        # Missing motion for a source outside the fitted stamp cannot invalidate
        # the in-stamp hypotheses. Included sources still fail closed.
        if not (-2<=x<cube.shape[2]+2 and -2<=y<cube.shape[1]+2): continue
        if status.startswith('missing'): astro=False
        fields=('ra_error','dec_error','pmra_error','pmdec_error')
        if any(n.get(v) is None for v in fields): astro=False
        else:
            epoch_year=2000+(float(np.median(t))+2454833-2451545)/365.25
            dt=epoch_year-float(n.get('ref_epoch',2016))
            sigma_mas=np.sqrt(sum(float(n[v])**2 for v in fields[:2])+dt**2*sum(float(n[v])**2 for v in fields[2:]))
            uncertainties[f'neighbour-{i}']=float(sigma_mas/(3980.))
            if uncertainties[f'neighbour-{i}']>.2: astro=False
        x,y=tpf.wcs.all_world2pix([[nra,ndec]],0)[0]
        if -2<=x<cube.shape[2]+2 and -2<=y<cube.shape[1]+2:
            positions[f'neighbour-{i}']=(origin[0]+x+.5,origin[1]+y+.5)
    if len(positions)<2: astro=False
    for name,pos in positions.items(): templates[name]=prf(cube.shape[1:],origin,pos)
    diff=difference_image(cube,t,candidate,tpf.pipeline_mask)
    selected_background=select_background_sources(neighbours)
    selection_path=out/'background-selection'/f"{hd.get('KEPLERID','unknown')}-q{quarter}.json"
    freeze(selection_path,{'seed':2026101303,'rules':'G 11-16; separation 2-30 arcsec; finite proper motion; no signal outcomes used',
        'sources':selected_background,'scope':'Development stamp unless independently selected pilot target'})
    provenance={'quarter':quarter,'module':module,'output':output,'prf_url':'https://archive.stsci.edu/missions/kepler/fpc/prf/'+prfname,
        'prf_sha256':source_hash(prfpath),'tpf_sha256':source_hash(path),'positions':positions,
        'gaia_position_uncertainty_pixels':uncertainties,'astrometric_uncertainty_method':'RA/Dec plus PM errors propagated in quadrature without covariances; >0.2 pixel unresolved',
        'calibration':'Official 2011265 channel calibration selected per quarter; no quarter-dependent focus calibration available'}
    if diff.get('status')!='screened': return {**provenance,'status':'unresolved','reason':diff.get('reason')}
    result=compare_hypotheses(diff['difference_flux'],diff['uncertainty'],templates,astrometry_complete=astro)
    masks={'pipeline':tpf.pipeline_mask,'all_pixels':np.ones(cube.shape[1:],bool)}
    aperture=[{'aperture':name,'difference_flux':float(np.nansum(diff['difference_flux'][mask])),
               'diagonal_error':float(np.sqrt(np.nansum(diff['uncertainty'][mask]**2)))} for name,mask in masks.items()]
    try:
        validation=extraction_control(templates) if len(templates)>=2 else {'passed':False,'reason':'No second hypothesis in stamp'}
    except ValueError as exc: validation={'passed':False,'reason':str(exc)}
    # Optional source extraction is enabled only after injection control validation.
    extracted=None
    if validation.get('passed'):
        f,e=extract_sources(cube,err,templates)
        csv_write(out/'extractions'/f"{path.stem}-q{quarter}.csv",pd.DataFrame({'time':t,**{name:f[:,j] for j,name in enumerate(templates)},
                    **{name+'_error':e[:,j] for j,name in enumerate(templates)}}))
        tf=pd.DataFrame({'time':t,'flux':f[:,0]/np.median(f[:,0]),'flux_err':e[:,0]/abs(np.median(f[:,0]))})
        tf=tf[np.isfinite(tf.flux)&np.isfinite(tf.flux_err)&(tf.flux_err>0)]
        known= fixed_holdout(tf,candidate) if candidate.get('benchmark_role') else {'passes_frozen_threshold':False,'reason':'No independently labelled source control'}
        validation['known_source_fixed_ephemeris']=known
        validated=bool(known.get('passes_frozen_threshold'))
        extracted={'status':'experimental_control_only','sources':list(templates),'known_source_timing_passed':validated,'real_background_discovery_use':False}
        # Both source injection and known-source checks are required by the optional API.
        if validated:
            extracted['validated_api_result_sources']=list(extract_background_discovery(cube,err,t,templates,list(templates)[1:],
                {'source_injection_passed':True,'known_source_timing_passed':True}))
    return serial({**provenance,**result,'aperture_comparison':aperture,'source_injection_validation':validation,'extraction':extracted})

# ---- pixel host test ----
MIN_FRAC = 0.25
QDATES = {1: '2009166043257', 2: '2009259160929', 3: '2009350155506', 4: '2010078095331', 5: '2010174085026',
          6: '2010265121752', 7: '2010355172524', 8: '2011073133259', 9: '2011177032512', 10: '2011271113734',
          11: '2012004120508', 12: '2012088054726'}
def write(p, v): p.parent.mkdir(parents=True, exist_ok=True); p.write_text(json.dumps(serial(v), indent=2, allow_nan=False))

class Downloads:
    def __init__(self, out):
        self.out = out
        self.records = read(out / 'downloads.json', [])
        self.total = sum(r['bytes'] for r in self.records)
    def download(self, url, path, max_bytes=100_000_000):
        path = pathlib.Path(path).resolve()
        if path.exists():  # offline reuse: a recorded file must still match its size and sha256
            rec = next((r for r in self.records if r['path'] == str(path.relative_to(self.out))), None)
            if rec and (path.stat().st_size != rec['bytes'] or source_hash(path) != rec['sha256']): raise ValueError(f'cached file changed: {path}')
            return path
        path.parent.mkdir(parents=True, exist_ok=True)
        part = path.with_name(path.name + '.partial')
        with urllib.request.urlopen(url, timeout=120) as r, part.open('wb') as s:
            data = r.read(max_bytes + 1)
            if len(data) > max_bytes or self.total + len(data) > 600_000_000: raise ValueError('download cap')
            s.write(data)
        part.replace(path); self.total += len(data)
        self.records.append({'url': url, 'path': str(path.relative_to(self.out)), 'bytes': len(data), 'sha256': source_hash(path)})
        write(self.out / 'downloads.json', self.records)
        return path

def read_leads(path):
    df = pd.read_csv(path, dtype={'gaia': str})
    missing = {'rank', 'lg', 'role', 'gaia', 'origin', 'peak_index', 'src'} - set(df.columns)
    if missing: raise ValueError(f'leads CSV missing columns: {sorted(missing)}')
    leads = df.to_dict(orient='records')
    for r in leads:
        for k in ('rank', 'lg', 'origin', 'peak_index'): r[k] = int(r[k])
        r['src'] = pathlib.Path(r['src'])
    return leads

def gaia_cone(lead, dest, limits):
    rc = read(lead['src'] / 'receipt.json')
    q = ("SELECT source_id,ra,dec,ref_epoch,pmra,pmdec,ra_error,dec_error,pmra_error,pmdec_error,phot_g_mean_mag "
         f"FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{rc['ra']},{rc['dec']},{35/3600}))")
    url = 'https://gea.esac.esa.int/tap-server/tap/sync?' + urllib.parse.urlencode({'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv', 'QUERY': q})
    f = pd.read_csv(limits.download(url, dest / 'gaia35.csv', 5_000_000), dtype={'source_id': str})
    host = f.loc[f.source_id == lead['gaia']].iloc[0]
    f['separation_arcsec'] = SkyCoord(f.ra.values * u.deg, f.dec.values * u.deg).separation(SkyCoord(host.ra * u.deg, host.dec * u.deg)).arcsec
    f['o'] = (f.source_id != lead['gaia']).astype(int)
    f = f.sort_values(['o', 'separation_arcsec']).drop(columns='o')
    return [{k: (None if pd.isna(v) else v) for k, v in r.items()} for r in f.to_dict(orient='records')]

def stamp_fractions(tpf, path, rows, prfdir, limits):
    """In-stamp PRF fraction per row, same geometry/PRF as localize_tpf."""
    hd = tpf.hdu[0].header; module, output = hd['MODULE'], hd['OUTPUT']
    name = f'kplr{module:02d}.{output}_2011265_prf.fits'
    prfpath = limits.download('https://archive.stsci.edu/missions/kepler/fpc/prf/' + name, prfdir / name, 30_000_000)
    shape = tpf.flux.shape[1:]
    prf = official_prf(prfpath, float(tpf.column) + shape[1] / 2, float(tpf.row) + shape[0] / 2, module, output)
    origin = (float(tpf.column), float(tpf.row)); tmid = float(np.median(tpf.time.bkjd))
    valid = np.isfinite(np.asarray(tpf.flux.value)).any(axis=0)
    out = []
    for r in rows:
        ra, dec, _ = propagate_gaia(r, tmid)
        x, y = tpf.wcs.all_world2pix([[ra, dec]], 0)[0]
        t = prf(shape, origin, (origin[0] + x + .5, origin[1] + y + .5))
        out.append({'x': float(x), 'y': float(y), 'frac': float(np.sum(t[valid]))})
    return out

def source_depths(csv, cand, names):
    f = pd.read_csv(csv); P, T0, D = cand['period_days'], cand['epoch_bkjd'], cand['duration_days']
    res = {}
    for n in names:
        flux = f[n].to_numpy(); t = f.time.to_numpy(); ok = np.isfinite(flux) & np.isfinite(t)
        def est(epoch, t=t, ok=ok, flux=flux):
            ph, ev, ins = transit_masks(t[ok], P, epoch, D); vals = []
            for e in np.unique(ev[ins]):
                i = ins & (ev == e); o = (ev == e) & (abs(ph) >= D / 2) & (abs(ph) < 2 * D)
                if i.sum() >= 1 and o.sum() >= 2:
                    b = np.mean(flux[ok][o]); vals.append((b - np.mean(flux[ok][i]), b))
            return vals
        res[n] = {'events': est(T0), 'null': [est(T0 + P * (k + .5) / 25) for k in range(24)]}
    return res

def host_A(per_q):
    assigned = [r['decision']['host'] for r in per_q if r.get('decision', {}).get('host')]
    if assigned:
        top = max(set(assigned), key=assigned.count)
        if assigned.count(top) >= 2 / 3 * len(assigned) and all(a == top for a in assigned): return top

def verdict_of(A, B, tgt):
    return 'KEEP' if A == tgt and B == tgt else 'DROP' if A and A == B and A != tgt else 'UNRESOLVED'

def process(lead, out, limits):
    import lightkurve as lk  # optional extra: pip install deimos[pixel]
    dest = out / str(lead['rank']); dest.mkdir(exist_ok=True)
    split = read(lead['src'] / 'split.json')['split']
    peak = read(lead['src'] / 'real.json')['peaks'][lead['peak_index']]
    cand = {'period_days': peak['period_days'], 'epoch_bkjd': peak['epoch_bkjd'], 'duration_days': peak['duration_days']}
    rows = gaia_cone(lead, dest, limits)
    quarters = [q for q in split['discovery_quarters'] if q in QDATES and q not in split['holdout_quarters']]
    k = f"{lead['origin']:09d}"; per_q = []
    for q in quarters:
        qfile = dest / f'q{q}.json'
        if qfile.exists(): per_q.append(read(qfile)); continue
        url = f'https://archive.stsci.edu/missions/kepler/target_pixel_files/{k[:4]}/{k}/kplr{k}-{QDATES[q]}_lpd-targ.fits.gz'
        try: path = limits.download(url, dest / 'tpf' / url.rsplit('/', 1)[1])
        except Exception as exc: per_q.append({'quarter': q, 'error': f'download: {exc}'}); continue
        tpf = lk.KeplerTargetPixelFile(path, quality_bitmask='default')
        assert int(tpf.quarter) == q and int(tpf.hdu[0].header['KEPLERID']) == lead['origin'] and q not in split['holdout_quarters']
        geo = stamp_fractions(tpf, path, rows, dest / 'prf', limits)
        keep = [r for r, g in zip(rows, geo) if r['separation_arcsec'] < 2 or g['frac'] >= MIN_FRAC]
        dropped = [{'source_id': r['source_id'], 'G': r['phot_g_mean_mag'], 'sep': r['separation_arcsec'], 'frac': g['frac']}
                   for r, g in zip(rows, geo) if r['separation_arcsec'] >= 2 and 0.005 < g['frac'] < MIN_FRAC]
        names = {'target': 'Gaia DR3 ' + keep[0]['source_id']}
        names.update({f'neighbour-{i}': 'Gaia DR3 ' + r['source_id'] for i, r in enumerate(keep)})
        res = localize_tpf(path, cand, keep, dest, limits)
        fits_ = res.get('fits', [])
        dec = {'host': None}
        if len(fits_) >= 2:
            d = fits_[1]['chi2'] - fits_[0]['chi2']
            ok = res['status'] in ('target_supported', 'off_target_supported') or d >= 9
            dec = {'host': names[fits_[0]['hypothesis']] if ok else None, 'best': names[fits_[0]['hypothesis']],
                   'second': names[fits_[1]['hypothesis']], 'delta_chi2': d, 'tool_status': res['status']}
        rec = {'quarter': q, 'decision': dec, 'names': names, 'dropped_edge_sources': dropped,
               'injection_passed': bool(res.get('source_injection_validation', {}).get('passed')),
               'injection': res.get('source_injection_validation'), 'localize_tpf': res}
        write(qfile, rec); per_q.append(rec); print(lead['rank'], q, dec, 'inj', rec['injection_passed'], flush=True)
    A = host_A(per_q)
    # Lead-level B: pool events across injection-passed quarters
    pooled = {}
    for r in per_q:
        if not r.get('injection_passed'): continue
        csvs = list((dest / 'extractions').glob(f"*{QDATES[r['quarter']]}*-q{r['quarter']}.csv"))
        if not csvs: continue
        cols = set(pd.read_csv(csvs[0], nrows=1).columns)
        for hyp, d in source_depths(csvs[0], cand, [n for n in r['names'] if n in cols]).items():
            nm = r['names'][hyp]; p = pooled.setdefault(nm, {'events': [], 'null': [[] for _ in range(24)], 'quarters': []})
            p['events'] += d['events']; p['quarters'].append(r['quarter'])
            for j in range(24): p['null'][j] += d['null'][j]
    B_table = []
    for nm, p in pooled.items():
        if not p['events']: continue
        frac = lambda ev: sum(a for a, b in ev) / sum(b for a, b in ev) if ev and sum(b for a, b in ev) > 0 else np.nan
        dep = frac(p['events']); nulls = np.array([frac(n) for n in p['null']]); nulls = nulls[np.isfinite(nulls)]
        sig = 1.4826 * np.median(abs(nulls - np.median(nulls))) if len(nulls) > 5 else np.nan
        meanflux = np.mean([b for a, b in p['events']])
        # A source with no positive baseline flux has an undefined depth: untestable, not a contradiction.
        z = float(dep / sig) if np.isfinite(dep) and np.isfinite(sig) and sig > 0 else None
        B_table.append({'source': nm, 'depth_ppm': dep * 1e6, 'null_sigma_ppm': sig * 1e6, 'z': z, 'measurable': z is not None,
                        'mean_flux_e_s': meanflux, 'dip_flux_e_s': dep * meanflux, 'n_events': len(p['events']), 'quarters': p['quarters']})
    B_table.sort(key=lambda r: -(r['z'] or -99))
    B = None
    strong = [r for r in B_table if (r['z'] or 0) >= 5]
    if len(strong) == 1 and all(r['z'] < 3 for r in B_table if r is not strong[0] and r['z'] is not None): B = strong[0]['source']
    tgt = 'Gaia DR3 ' + lead['gaia']
    verdict = verdict_of(A, B, tgt)
    s = {'rank': lead['rank'], 'role': lead['role'], 'target': tgt, 'candidate': cand, 'quarters_used': quarters,
         'per_quarter': [{'quarter': r['quarter'], **r.get('decision', {}), 'injection_passed': r.get('injection_passed'),
                          'error': r.get('error'), 'n_dropped_edge': len(r.get('dropped_edge_sources', []))} for r in per_q],
         'A_difference_image_host': A, 'B_extraction_host': B, 'B_table': B_table, 'verdict': verdict}
    write(dest / 'summary.json', s)
    return s

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('leads_csv'); ap.add_argument('--out', required=True); ap.add_argument('ranks', nargs='*', type=int)
    a = ap.parse_args(argv)
    out = pathlib.Path(a.out).resolve()
    for name in ['cache', 'tmp', 'config']: (out / name).mkdir(parents=True, exist_ok=True)
    os.environ.update({'XDG_CACHE_HOME': str(out / 'cache'), 'XDG_CONFIG_HOME': str(out / 'config'), 'MPLCONFIGDIR': str(out / 'config'),
                       'TMPDIR': str(out / 'tmp'), 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1'})
    leads = read_leads(a.leads_csv)
    limits = Downloads(out)
    allres = []
    for lead in leads:
        if a.ranks and lead['rank'] not in a.ranks: continue
        try: allres.append(process(lead, out, limits))
        except Exception:
            (out / f"error-{lead['rank']}.txt").write_text(traceback.format_exc()); print(traceback.format_exc(), flush=True)
    print(json.dumps(serial([{k: v for k, v in s.items() if k in ('rank', 'A_difference_image_host', 'B_extraction_host', 'verdict')} for s in allres]), indent=1))

if __name__ == '__main__':
    main()
