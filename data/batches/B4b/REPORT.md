# B4b: paired Kepler Bonus search

Status: STOPPED_USER. Stop: User ended cloud search; preserve completed work; never replay interrupted attempts. 2374/5260 sources terminal; 2287 complete, 87 failed, 0 interrupted (never replayed).

```json
{
  "real": {
    "fraction": 1.0,
    "completed_arms": 2287,
    "survivor_peaks": 1634,
    "survivor_sources": 687,
    "survivor_dependency_groups": 41,
    "cache_hits": 0
  },
  "inverted": {
    "fraction": 1.0,
    "completed_arms": 2287,
    "survivor_peaks": 1686,
    "survivor_sources": 713,
    "survivor_dependency_groups": 41,
    "cache_hits": 0
  },
  "block_scrambled": {
    "fraction": 1.0,
    "completed_arms": 2287,
    "survivor_peaks": 427,
    "survivor_sources": 243,
    "survivor_dependency_groups": 36,
    "cache_hits": 0
  },
  "injected": {
    "fraction": 1.0,
    "completed_arms": 2287,
    "survivor_peaks": 3388,
    "survivor_sources": 1319,
    "survivor_dependency_groups": 41,
    "cache_hits": 0
  }
}
```

Injection detected/recovered: 1340/890 of 2287 completed sources. Lead groups: 1495; holdout-supported groups: 0; unresolved: 1495.

Source-level denominators are correlated; dependency groups are conservative and not certified independent. Real/control contrasts describe workload, not discovery yield. One flux-level injection per source cannot validate extraction. Host assignment, completeness, FPP and novelty remain unresolved; no promotion. Rows: reports/*.jsonl (one row per source per arm; peaks.jsonl one row per peak).
