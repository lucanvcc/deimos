"""Per-source arm science: input transforms, injection, recovery.

Everything in this file (plus the adapter, the search plug-in and core.py) is covered by
code_hash(), which is recorded in every receipt.
"""
import sys, pathlib
import numpy as np
from .core import physical_transit
from .core import transit_masks
from .util import digest, sha


def transform(d, kind, seed):
 out=d.copy().reset_index(drop=True)
 if kind=='inverted':out['flux']=2-out.flux
 elif kind=='block_scrambled':
  rng=np.random.default_rng(seed)
  for _,indices in out.groupby('quarter').groups.items():
   indices=np.asarray(indices);tt=out.loc[indices,'time'].to_numpy();cuts=np.r_[0,np.flatnonzero(np.diff(tt)>.5)+1,len(tt)]
   for a,b in zip(cuts[:-1],cuts[1:]):
    ix=indices[a:b]
    if len(ix)<3:continue
    cadence=np.median(np.diff(tt[a:b]));size=max(2,int(round(1/cadence)));n=len(ix)//size
    if n<3:continue
    blocks=ix[:n*size].reshape(n,size);out.loc[blocks.ravel(),['flux','flux_err']]=out.loc[blocks[rng.permutation(n)].ravel(),['flux','flux_err']].to_numpy()
 elif kind!='real':raise ValueError(kind)
 return out


def truth_for(rank, d, seed_base, injection):
 """Injection truth for one source; the grid comes from the manifest."""
 rng=np.random.default_rng(seed_base+rank);lo,hi=injection['jitter']
 nP,nD,nW=len(injection['periods_days']),len(injection['depths']),len(injection['durations_days'])
 idx=(rank-1)%(nP*nD*nW);cell=idx
 p=injection['periods_days'][idx//(nD*nW)]*rng.uniform(lo,hi);depth=injection['depths'][(idx//nW)%nD]*rng.uniform(lo,hi);duration=injection['durations_days'][idx%nW]*rng.uniform(lo,hi)
 return dict(period_days=p,depth=depth,duration_days=duration,epoch_bkjd=float(d.time.min()+rng.uniform(0,p)),cell=cell)


def matched(peak,truth,times):
 _,_,a=transit_masks(times,peak['period_days'],peak['epoch_bkjd'],peak['duration_days']);_,_,b=transit_masks(times,truth['period_days'],truth['epoch_bkjd'],truth['duration_days'])
 return abs(peak['period_days']/truth['period_days']-1)<.01 and a.any() and b.any() and (a&b).sum()/max(1,b.sum())>=.5


def recovered(p,t,times):return bool(matched(p,t,times) and 1/1.5<=p['duration_days']/t['duration_days']<=1.5)


def arm_seed(seed_base, rank):
    return seed_base + int(rank)


def arm_frame(d, arm, rank, seed_base, truth, injection):
    """The searched frame for one arm. """
    x = d.copy() if arm == 'injected' else transform(d, arm, arm_seed(seed_base, rank))
    if arm == 'injected':
        x['flux'] *= physical_transit(x.time.to_numpy(), truth['period_days'], truth['epoch_bkjd'], truth['duration_days'], truth['depth'], exposure=injection['exposure_days'])
    return x


def code_hash(plugin_module):
    """Hash of all code that can change an arm result: science.py, adapter, plug-in, core.py, library versions."""
    here = pathlib.Path(__file__).resolve().parent
    files = [here / 'science.py', here / 'adapter.py', here / 'core.py', pathlib.Path(plugin_module.__file__).resolve()]
    import numpy, pandas, astropy
    return digest(dict(files={str(p.name if p.parent == here else p.relative_to(p.parents[1])): sha(p) for p in files},
                       python=sys.version.split()[0], numpy=numpy.__version__, pandas=pandas.__version__,
                       astropy=astropy.__version__))
