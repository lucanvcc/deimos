"""Numerical search plug-ins.

A plug-in is any callable ``select(frame, cfg) -> dict`` returning the arm packet
(peaks with diagnostics, harmonic_groups, seconds, rows, baseline_days). Name it in
the manifest as "module:function". Its source file is part of the code hash, so
swapping it invalidates every cached arm result automatically.
"""
import importlib, time
import numpy as np
from .core import clean_frame
from .core import extract_multi
from .core import diagnostics, serial


def frozen_select(raw, cfg):
    """Standard search: clean, multi-peak BLS, diagnostics."""
    t = time.monotonic(); clean = clean_frame(raw, cfg['clean']); multi = extract_multi(clean, cfg['search'], max_peaks=5, residual=True); peaks = list(multi['peaks'])
    if multi['residual'] and multi['residual'].get('status') == 'ok': peaks.append(multi['residual'])
    for p in peaks: p['diagnostics'] = diagnostics(clean, p, cfg['validate'])
    return serial(dict(peaks=peaks, harmonic_groups=multi['harmonic_groups'], seconds=time.monotonic() - t, rows=len(clean), baseline_days=float(np.ptp(clean.time))))


def load_plugin(spec):
    module, _, name = spec.partition(':')
    mod = importlib.import_module(module)
    return getattr(mod, name), mod
