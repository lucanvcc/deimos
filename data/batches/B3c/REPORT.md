# B3c: paired Kepler Bonus search

Status: COMPLETE. Stop: None. 5000/5000 sources terminal; 4900 complete, 100 failed, 0 interrupted (never replayed).

```json
{
  "real": {
    "fraction": 1.0,
    "completed_arms": 4900,
    "survivor_peaks": 3748,
    "survivor_sources": 1615,
    "survivor_dependency_groups": 42,
    "cache_hits": 0
  },
  "inverted": {
    "fraction": 1.0,
    "completed_arms": 4900,
    "survivor_peaks": 3877,
    "survivor_sources": 1707,
    "survivor_dependency_groups": 43,
    "cache_hits": 0
  },
  "block_scrambled": {
    "fraction": 1.0,
    "completed_arms": 4900,
    "survivor_peaks": 897,
    "survivor_sources": 503,
    "survivor_dependency_groups": 37,
    "cache_hits": 0
  },
  "injected": {
    "fraction": 1.0,
    "completed_arms": 4900,
    "survivor_peaks": 7213,
    "survivor_sources": 2835,
    "survivor_dependency_groups": 44,
    "cache_hits": 0
  }
}
```

Injection detected/recovered: 2769/1811 of 4900 completed sources. Lead groups: 3312; holdout-supported groups: 80; unresolved: 362.

Source-level denominators are correlated; dependency groups are conservative and not certified independent. Real/control contrasts describe workload, not discovery yield. One flux-level injection per source cannot validate extraction. Host assignment, completeness, FPP and novelty remain unresolved; no promotion. Rows: reports/*.jsonl (one row per source per arm; peaks.jsonl one row per peak).
