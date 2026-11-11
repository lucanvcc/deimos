"""Fast deterministic stand-in for the numerical search, used only by tests."""
import time
import numpy as np


def select(frame, cfg):
    t = time.monotonic()
    flux = frame.flux.to_numpy()
    snr = float(abs(np.mean(flux) - 1) / (np.std(flux) / np.sqrt(len(flux))) + 8 * (flux.argmin() % 2))
    peak = dict(period_days=2.0 + float(flux.argmin() % 7), epoch_bkjd=float(frame.time.min()) + 0.3, duration_days=0.1, depth=float(1 - flux.min()),
                snr=snr, number_of_transits=5, initial_candidate=snr > 7, rank=1,
                diagnostics=dict(photometry_vetting_pass=bool(snr > 7)))
    return dict(peaks=[peak], harmonic_groups=[], seconds=time.monotonic() - t, rows=len(frame), baseline_days=float(np.ptp(frame.time)))
