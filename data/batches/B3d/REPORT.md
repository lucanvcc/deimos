# B3d: paired Kepler Bonus search

Status: COMPLETE. Stop: None. 5000/5000 sources terminal; 4776 complete, 224 failed, 0 interrupted (never replayed).

```json
{
  "real": {
    "fraction": 1.0,
    "completed_arms": 4776,
    "survivor_peaks": 3192,
    "survivor_sources": 1373,
    "survivor_dependency_groups": 42,
    "cache_hits": 0
  },
  "inverted": {
    "fraction": 1.0,
    "completed_arms": 4776,
    "survivor_peaks": 3394,
    "survivor_sources": 1469,
    "survivor_dependency_groups": 42,
    "cache_hits": 0
  },
  "block_scrambled": {
    "fraction": 1.0,
    "completed_arms": 4776,
    "survivor_peaks": 877,
    "survivor_sources": 466,
    "survivor_dependency_groups": 37,
    "cache_hits": 0
  },
  "injected": {
    "fraction": 1.0,
    "completed_arms": 4776,
    "survivor_peaks": 6770,
    "survivor_sources": 2676,
    "survivor_dependency_groups": 43,
    "cache_hits": 0
  }
}
```

Injection detected/recovered: 2686/1778 of 4776 completed sources. Lead groups: 2805; holdout-supported groups: 91; unresolved: 371.

Source-level denominators are correlated; dependency groups are conservative and not certified independent. Real/control contrasts describe workload, not discovery yield. One flux-level injection per source cannot validate extraction. Host assignment, completeness, FPP and novelty remain unresolved; no promotion. Rows: reports/*.jsonl (one row per source per arm; peaks.jsonl one row per peak).
