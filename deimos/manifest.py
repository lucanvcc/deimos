"""Manifest loading and validation. A batch is fully described by its manifest.

Unknown fields are rejected so that a typo cannot silently change what is run.
"""
import os, pathlib, json
import pandas as pd
from .util import digest, load, sha, ENV_STATE, ENV_INPUT

KNOWN_ARMS = ('real', 'inverted', 'block_scrambled', 'injected')

DEFAULT_INJECTION = {
    # cell=(rank-1)%8 selects period/depth/duration, each jittered by U(.9,1.1).
    'periods_days': [2.0, 10.0],
    'depths': [0.0005, 0.002],
    'durations_days': [0.1, 0.2],
    'jitter': [0.9, 1.1],
    'exposure_days': 1765.464 / 86400,
}

DEFAULT_SCIENCE_CONFIG = pathlib.Path(__file__).resolve().parent / 'science-config.json'
DEFAULT_PLUGIN = 'deimos.search:frozen_select'

REQUIRED = ['registration_csv', 'output_dir', 'seed_base']
SCIENTIFIC = ['name', 'registration_csv', 'science_config', 'seed_base', 'select', 'stages', 'holdout_stage', 'holdouts',
              'arms', 'injection', 'search_plugin', 'ranking_rule', 'scientific_method_signoff']
OPERATIONAL = ['output_dir', 'workers', 'threads_per_worker', 'export_every_seconds', 'state_dir', 'input_dirs']
KNOWN_FIELDS = set(SCIENTIFIC + OPERATIONAL)
REGISTRATION_COLUMNS = ['rank', 'gaia_designation', 'conservative_group', 'tpforg', 'kic', 'ra', 'dec', 'file_name', 'url']


class ManifestError(ValueError):
    """A manifest or registration file that cannot be run as written."""


def _need(cond, message):
    if not cond:
        raise ManifestError(message)


class Batch:
    """Validated, path-resolved view of one manifest."""

    def __init__(self, path):
        self.path = pathlib.Path(path).resolve()
        _need(self.path.is_file(), f'manifest not found: {path}')
        try:
            raw = json.loads(self.path.read_text())
        except json.JSONDecodeError as e:
            raise ManifestError(f'manifest is not valid JSON: {e}') from None
        _need(isinstance(raw, dict), 'manifest must be a JSON object')
        unknown = sorted(set(raw) - KNOWN_FIELDS)
        _need(not unknown, f'unknown manifest fields {unknown}; known: {sorted(KNOWN_FIELDS)}')
        missing = [k for k in REQUIRED if k not in raw]
        _need(not missing, f'batch manifest missing fields: {missing}')
        base = self.path.parent
        resolve = lambda p: (base / p).resolve() if p else None
        self.raw = raw
        self.name = str(raw.get('name', self.path.stem))
        self.registration_csv = resolve(raw['registration_csv'])
        self.output_dir = resolve(raw['output_dir'])
        self.science_config_path = resolve(raw['science_config']) if raw.get('science_config') else DEFAULT_SCIENCE_CONFIG
        _need(isinstance(raw['seed_base'], int) and not isinstance(raw['seed_base'], bool), 'seed_base must be an integer')
        self.seed_base = int(raw['seed_base'])
        self.select_spec = raw.get('select', {})
        _need(isinstance(self.select_spec, dict), 'select must be an object')
        bad = sorted(set(self.select_spec) - {'ranks', 'min_rank', 'max_rank', 'exclude_ranks', 'role'})
        _need(not bad, f'unknown select keys {bad}')
        self.holdouts = bool(raw.get('holdouts', False))
        try:
            self.workers = int(raw.get('workers', 1))
            self.threads_per_worker = int(raw.get('threads_per_worker', 1))
            self.export_every_seconds = float(raw.get('export_every_seconds', 120))
        except (TypeError, ValueError) as e:
            raise ManifestError(f'operational fields must be numbers: {e}') from None
        _need(self.workers >= 1, 'workers must be >= 1')
        self.injection = {**DEFAULT_INJECTION, **raw.get('injection', {})}
        self.search_plugin = raw.get('search_plugin', DEFAULT_PLUGIN)
        _need(isinstance(self.search_plugin, str) and ':' in self.search_plugin, 'search_plugin must look like "module:function"')
        self.arms = []
        _need(isinstance(raw.get('arms', list(KNOWN_ARMS)), list) and raw.get('arms', list(KNOWN_ARMS)), 'arms must be a non-empty list')
        for arm in raw.get('arms', list(KNOWN_ARMS)):
            arm = {'name': arm} if isinstance(arm, str) else dict(arm)
            _need(arm.get('name') in KNOWN_ARMS, f"unknown arm {arm.get('name')!r}; known: {KNOWN_ARMS}")
            _need(set(arm) == {'name'}, f"arm {arm['name']} takes only a name; per-arm options were removed")
            self.arms.append(arm)
        names = [a['name'] for a in self.arms]
        _need(len(set(names)) == len(names), 'duplicate arm names')
        _need(not self.holdouts or 'real' in names, 'holdouts need the real arm (leads come from real)')
        # Operational locations: environment (set by CLI flags) > manifest > default under the output dir.
        self.state_dir = pathlib.Path(os.environ[ENV_STATE]).resolve() if os.environ.get(ENV_STATE) else resolve(raw.get('state_dir'))
        _need(not self.holdouts or self.state_dir is not None,
              'holdouts: true needs an explicit state_dir (--state-dir, DEIMOS_STATE or "state_dir"): the holdout ledger must live '
              'outside this batch folder so a second batch cannot reopen the same sealed holdout')
        if self.state_dir is None:
            self.state_dir = self.output_dir / 'state'
        env_in = [p for p in os.environ.get(ENV_INPUT, '').split(os.pathsep) if p]
        self.input_dirs = [pathlib.Path(p).resolve() for p in env_in] + [resolve(p) for p in raw.get('input_dirs', [])]
        if 'stages' in raw:
            self.stages = [int(s) for s in raw['stages']]
            _need(self.stages == sorted(set(self.stages)) and self.stages and self.stages[0] > 0, 'stages must be strictly increasing positive source counts')
        else:
            self.stages = [None]     # filled by validate() once the registration is read
        self.holdout_stage_raw = raw.get('holdout_stage')

    # -- validation that needs the files -----------------------------------------------------------
    def validate(self):
        """Read and check the registration CSV, science config and selection. Returns the selected rows."""
        _need(self.registration_csv.is_file(), f'registration_csv not found: {self.registration_csv.name}')
        _need(self.science_config_path.is_file(), f'science_config not found: {self.science_config_path.name}')
        cfg = self.science_config()
        for section in ('clean', 'search', 'validate'):
            _need(isinstance(cfg.get(section), dict), f'science_config is missing the "{section}" section')
        reg = self.registration()
        missing = [c for c in REGISTRATION_COLUMNS if c not in reg.columns]
        _need(not missing, f'registration CSV missing columns: {missing}')
        _need(len(reg) > 0, 'registration CSV has no rows')
        try:
            ranks = reg['rank'].astype(int)
        except (TypeError, ValueError):
            raise ManifestError('registration ranks must be integers') from None
        _need(ranks.is_unique and (ranks > 0).all(), 'registration ranks must be unique positive integers')
        for col in ('ra', 'dec'):
            _need(pd.to_numeric(reg[col], errors='coerce').notna().all(), f'registration column {col} must be numeric')
        _need(reg['ra'].between(0, 360).all() and reg['dec'].between(-90, 90).all(), 'registration ra/dec out of range')
        for name in reg['file_name']:
            _need(isinstance(name, str) and name == pathlib.Path(name).name and name.endswith('.fits'), f'file_name must be a bare *.fits name, got {name!r}')
        if 'sha256' in reg.columns:
            bad = [s for s in reg['sha256'].dropna() if len(str(s)) != 64 or set(str(s)) - set('0123456789abcdef')]
            _need(not bad, f'registration sha256 values must be 64 lowercase hex characters: {bad[:2]}')
        rows = self._select(reg)
        if self.stages == [None]:
            self.stages = [len(rows)]
        _need(len(rows) == self.stages[-1], f'selection has {len(rows)} sources but final stage is {self.stages[-1]}; adjust "select" or "stages" so they agree')
        self.holdout_stage = int(self.holdout_stage_raw) if self.holdout_stage_raw is not None else self.stages[-1]
        _need(self.holdout_stage in self.stages, 'holdout_stage must be one of stages')
        return rows

    def science_config(self):
        return load(self.science_config_path)

    def config_hash(self):
        """Hash of every input that changes numerical arm output, excluding code (see science.code_hash)."""
        cfg = self.science_config()
        return digest(dict(clean=cfg['clean'], search=cfg['search'], validate=cfg['validate'],
                           injection=self.injection, search_plugin=self.search_plugin,
                           search_strategy={'name': 'exact'}))

    def registration(self):
        return pd.read_csv(self.registration_csv)

    def _select(self, frame):
        frame = frame.sort_values('rank', kind='stable')
        spec = self.select_spec
        if 'ranks' in spec:
            wanted = [int(r) for r in spec['ranks']]
            unknown = sorted(set(wanted) - set(frame['rank'].astype(int)))
            _need(not unknown, f'select.ranks not in the registration: {unknown}')
            frame = frame.set_index('rank', drop=False).loc[wanted].reset_index(drop=True)
        if 'min_rank' in spec:
            frame = frame[frame['rank'] >= int(spec['min_rank'])]
        if 'max_rank' in spec:
            frame = frame[frame['rank'] <= int(spec['max_rank'])]
        if spec.get('exclude_ranks'):
            frame = frame[~frame['rank'].isin([int(r) for r in spec['exclude_ranks']])]
        if spec.get('role'):
            _need('role' in frame.columns, 'select.role needs a "role" column in the registration CSV')
            frame = frame[frame['role'].isin(spec['role'] if isinstance(spec['role'], list) else [spec['role']])]
        return frame.to_dict('records')

    def selection(self):
        """Ordered rows to search (validates the registration on every call; it is small)."""
        return self.validate()

    def identity(self):
        """Scientific identity written into the batch folder; a resume must match it exactly.
        Operational knobs (workers, export interval, state/input locations) may change between resumes."""
        self.selection()
        return dict(registration_sha256=sha(self.registration_csv), science_config_sha256=sha(self.science_config_path),
                    config_hash=self.config_hash(), seed_base=self.seed_base, stages=self.stages,
                    holdout_stage=self.holdout_stage, holdouts=self.holdouts, arms=self.arms, select=self.select_spec)
