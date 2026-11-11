# Method

This describes what the code in `deimos/` does. It is the procedure that produced the candidates in
`data/candidates/`. It finds and vets candidate signals; confirming a planet takes follow-up outside this repo.

## 1. Inputs and adapter (`adapter.py`)

- Input is a Kepler Bonus `kbonus-bkg` FITS file (primary header `HLSPID` must be `KBONUS-BKG`), found on local disk by
  `file_name`. If the registration has a `sha256` column, a file with a different digest is refused.
- The Gaia designation in the primary header (`GAIAID`) must match the registration row.
- Each `LIGHTCURVE_Q<n>` extension is read (`TIME`, `FLUX`, `FLUX_ERR`). `FLUX` and `FLUX_ERR` must be in `e-/s`,
  `TIME` in `BJD - 2454833` with `TIMESYS` TDB. Anything else raises an error.
- Cadences are kept only if time, flux and error are finite, the error is positive, and `SAP_QUALITY == 0`. A quarter
  whose header `PSFFRAC` is not above 0.5 is dropped.
- Flux and error are divided by the quarter median. Duplicate times are removed. Time is BKJD (BJD_TDB - 2454833).
- The configured flux choice and quality policy are recorded in `science-config.json` (`adapter_flux`, `quality_mask`).
  The quality policy is deliberately strict and is not assumed equivalent to the archive's default masks.

## 2. Discovery / holdout quarter split (`core.py: split_quarters`)

- Quarters with at least 1000 usable cadences are eligible; at least 4 are required or the source fails.
- The last `ceil(25%)` of eligible quarters are reserved as holdout; the rest are the discovery set. For the
  example sources this gives 7 discovery and 3 holdout quarters.
- The split uses coverage only, before any flux is searched. Search arms only ever see discovery quarters.
- The holdout quarters are written to `sealed-holdout.csv.gz` only if the manifest has `holdouts: true`. With
  `holdouts: false` they are not written and not used.
- `upstream_independence` is recorded as false: the holdout quarters share upstream processing (a PSF model fit to
  quarter-wide pixels and a shared pixel stack) with the discovery quarters.

## 3. Arms (`science.py`)

All arms use the same discovery frame and the same search. `arm_seed = seed_base + rank`.

| arm | construction |
|---|---|
| `real` | the discovery light curve unchanged |
| `inverted` | flux replaced by `2 - flux` (reflected about the normalised level of 1) |
| `block_scrambled` | within each quarter and each gap-free segment (gaps > 0.5 d), blocks of about 1 day (rounded to whole cadences) are shuffled with `numpy.random.default_rng(seed)`; flux and error move together; short segments are left unchanged |
| `injected` | a simulated transit is multiplied into the real light curve (`physical_transit`: opaque circular body on a uniform stellar disk with no limb darkening, impact parameter 0.3, duration taken as first-to-fourth contact, 21-point supersampling over the exposure time) |

Injection parameters per source come from `rng = default_rng(seed_base + rank)`: the grid cell is
`(rank - 1) mod (n_periods * n_depths * n_durations)`, and period, depth and duration are each multiplied by a draw
from U(0.9, 1.1); the epoch is drawn uniformly within one period of the first discovery time. The truth is written
to `private-injection-truth.json` in the source folder. The injected signal is synthetic by construction.

Recovery (`injection_detected`, `injection_recovered`): a peak matches the injection if its period is within 1 %,
at least half of the injected in-transit samples fall inside the peak's transit mask, and the duration ratio is
within 1/1.5 to 1.5. "Detected" additionally needs `initial_candidate`; "recovered" needs `photometry_vetting_pass`.

## 4. Frozen clean, search and vetting (`search.py`, `core.py`)

Settings are in `science-config.json`.

- **Clean** (`clean_frame`): per quarter and gap-free segment (gap 0.5 d), divide by the median, flatten with a
  running median of 2.1 d, clip at +5 sigma and -20 sigma (robust scatter), require at least 100 cadences.
- **Search** (`extract_bls`, `extract_multi`): astropy `BoxLeastSquares`, SNR objective, periods 0.5 to 30 d on a
  geometric grid whose spacing limits phase drift to half the shortest duration, 8 durations from 0.05 to 0.3 d.
  Errors are floored by robust point-to-point scatter. Up to 5 distinct local maxima are retained (peaks within 1 %
  of 0.25, 0.5, 1, 2 or 4 times a kept period are grouped as harmonics), plus one more search after masking the
  strongest peak at 1.5 times its duration.
- **Initial candidate**: SNR at least 7, at least 3 observed transits with at least 2 points each, depth between 0
  and 1. SNR is a white-noise depth/error statistic and is not corrected for the number of trials.
- **Vetting** (`diagnostics`): odd/even depth difference at most 3 sigma, secondary below 5 sigma, no single event
  above 80 % of the signal, remaining-event SNR at least 3, plus a repair screen (red-noise corrected SNR at least 7,
  no broad modulation, at least 3 measurable events). The result is `photometry_vetting_pass`. Localisation is
  recorded as `unresolved`; no pixel-level source assignment is done.

A "survivor" in the outputs is a peak with `photometry_vetting_pass` true.

## 5. Per-source receipts, stages and lead groups

- Each source runs in its own worker process. `attempt.json` is written first, `receipt.json` last. A source with an
  attempt but no receipt (interrupted) is marked `interrupted_incomplete` on resume and never replayed.
- Sources are ordered by registration `rank`. A **stage** is a prefix of that order (`stages` in the manifest;
  default one stage with all sources). When every source in a prefix has a receipt, `stage-N-freeze.json` freezes the
  ranked real-arm survivors: sorted by descending SNR, then source designation, then peak index.
- **Lead groups** (`report.lead_groups`): real-arm survivors are merged within the same `conservative_group` when
  they share an ephemeris (period ratio within 1 % of 0.25, 0.5, 1, 2 or 4 and at least half of the sampled transit
  points overlapping). `host_status` and `novelty_status` are `unresolved` for every lead.

## 6. Holdout check (optional; `stages.py`)

Only with `holdouts: true`, which requires an explicit state directory outside the batch folder.

- Evaluated once, after the holdout stage is frozen, only for leads in the frozen list, at the frozen period, epoch
  and duration. No refitting. The measurement is a fixed-template depth with a 7-sigma and 3-transit threshold.
- An access record (`access-<peak>.json`) is written before the sealed file is opened. An access without a result is
  never retried (`interrupted_holdout_unresolved`).
- A ledger in `<state_dir>/holdout-ledger/`, keyed by the input FITS sha256, refuses a second batch that wants the
  same source (`holdout_previously_accessed_by_other_batch`).
- Statuses: `supported_fixed_photometry_host_unresolved`, `failed_frozen_screen`, `insufficient_power`,
  `interrupted_holdout_unresolved`, `holdout_previously_accessed_by_other_batch`. Even a supported lead has an
  unresolved host and is not validated.
- Holdout independence from upstream processing is not established.

## 7. Frozen pins and reproducibility

- `frozen-pins.json` holds SHA-256 digests of `core.py` and `adapter.py`. `run`, workers and `verify-pins`
  fail if any file differs.
- `receipt.json` and `provenance.json` record `config_hash` (clean/search/validate settings, injection grid, plug-in),
  `code_hash` (science, adapter, search plug-in, core.py, Python/numpy/pandas/astropy versions), seeds, input
  sha256 and package version.
- A resume must match `batch-identity.json` (registration, science config, seed, stages, arms, selection) or it is
  refused; use a new `output_dir` for a different configuration.

## 8. Not implemented

- Anything other than Kepler Bonus `kbonus-bkg` FITS; TESS; Kepler PDCSAP or target pixel files.
- Host or pixel-level localisation, centroid or difference-image checks, catalogue cross-matches, false-positive
  probability, completeness, novelty assessment.
- Non-exact search strategies or alternative fast backends (a manifest asking for them is rejected).
- Downloads inside `run` (only `fetch` downloads).
- Windows support.
