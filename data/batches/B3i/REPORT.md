# B3i: paired Kepler Bonus search

Status: COMPLETE. Stop: None. 5000/5000 sources terminal; 4096 complete, 904 failed, 0 interrupted (never replayed).

```json
{
  "real": {
    "fraction": 1.0,
    "completed_arms": 4096,
    "survivor_peaks": 3348,
    "survivor_sources": 1445,
    "survivor_dependency_groups": 47,
    "cache_hits": 0
  },
  "inverted": {
    "fraction": 1.0,
    "completed_arms": 4096,
    "survivor_peaks": 3605,
    "survivor_sources": 1554,
    "survivor_dependency_groups": 43,
    "cache_hits": 0
  },
  "block_scrambled": {
    "fraction": 1.0,
    "completed_arms": 4096,
    "survivor_peaks": 827,
    "survivor_sources": 425,
    "survivor_dependency_groups": 38,
    "cache_hits": 0
  },
  "injected": {
    "fraction": 1.0,
    "completed_arms": 4096,
    "survivor_peaks": 5609,
    "survivor_sources": 2264,
    "survivor_dependency_groups": 48,
    "cache_hits": 0
  }
}
```

Injection detected/recovered: 1950/1253 of 4096 completed sources. Lead groups: 3019; holdout-supported groups: 55; unresolved: 411.

Source-level denominators are correlated; dependency groups are conservative and not certified independent. Real/control contrasts describe workload, not discovery yield. One flux-level injection per source cannot validate extraction. Host assignment, completeness, FPP and novelty remain unresolved; no promotion. Rows: reports/*.jsonl (one row per source per arm; peaks.jsonl one row per peak).
