# B3g: paired Kepler Bonus search

Status: COMPLETE. Stop: None. 5000/5000 sources terminal; 4603 complete, 397 failed, 0 interrupted (never replayed).

```json
{
  "real": {
    "fraction": 1.0,
    "completed_arms": 4603,
    "survivor_peaks": 3812,
    "survivor_sources": 1633,
    "survivor_dependency_groups": 52,
    "cache_hits": 0
  },
  "inverted": {
    "fraction": 1.0,
    "completed_arms": 4603,
    "survivor_peaks": 3997,
    "survivor_sources": 1733,
    "survivor_dependency_groups": 52,
    "cache_hits": 0
  },
  "block_scrambled": {
    "fraction": 1.0,
    "completed_arms": 4603,
    "survivor_peaks": 855,
    "survivor_sources": 481,
    "survivor_dependency_groups": 39,
    "cache_hits": 0
  },
  "injected": {
    "fraction": 1.0,
    "completed_arms": 4603,
    "survivor_peaks": 6684,
    "survivor_sources": 2643,
    "survivor_dependency_groups": 53,
    "cache_hits": 0
  }
}
```

Injection detected/recovered: 2392/1561 of 4603 completed sources. Lead groups: 3335; holdout-supported groups: 93; unresolved: 384.

Source-level denominators are correlated; dependency groups are conservative and not certified independent. Real/control contrasts describe workload, not discovery yield. One flux-level injection per source cannot validate extraction. Host assignment, completeness, FPP and novelty remain unresolved; no promotion. Rows: reports/*.jsonl (one row per source per arm; peaks.jsonl one row per peak).
