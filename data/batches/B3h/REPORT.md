# B3h: paired Kepler Bonus search

Status: COMPLETE. Stop: None. 5000/5000 sources terminal; 4326 complete, 674 failed, 0 interrupted (never replayed).

```json
{
  "real": {
    "fraction": 1.0,
    "completed_arms": 4326,
    "survivor_peaks": 3711,
    "survivor_sources": 1605,
    "survivor_dependency_groups": 43,
    "cache_hits": 0
  },
  "inverted": {
    "fraction": 1.0,
    "completed_arms": 4326,
    "survivor_peaks": 3922,
    "survivor_sources": 1659,
    "survivor_dependency_groups": 43,
    "cache_hits": 0
  },
  "block_scrambled": {
    "fraction": 1.0,
    "completed_arms": 4326,
    "survivor_peaks": 854,
    "survivor_sources": 435,
    "survivor_dependency_groups": 40,
    "cache_hits": 0
  },
  "injected": {
    "fraction": 1.0,
    "completed_arms": 4326,
    "survivor_peaks": 6019,
    "survivor_sources": 2430,
    "survivor_dependency_groups": 44,
    "cache_hits": 0
  }
}
```

Injection detected/recovered: 2052/1332 of 4326 completed sources. Lead groups: 3328; holdout-supported groups: 102; unresolved: 452.

Source-level denominators are correlated; dependency groups are conservative and not certified independent. Real/control contrasts describe workload, not discovery yield. One flux-level injection per source cannot validate extraction. Host assignment, completeness, FPP and novelty remain unresolved; no promotion. Rows: reports/*.jsonl (one row per source per arm; peaks.jsonl one row per peak).
