# B3j: paired Kepler Bonus search

Status: COMPLETE. Stop: None. 5000/5000 sources terminal; 3860 complete, 1140 failed, 0 interrupted (never replayed).

```json
{
  "real": {
    "fraction": 1.0,
    "completed_arms": 3860,
    "survivor_peaks": 3137,
    "survivor_sources": 1352,
    "survivor_dependency_groups": 45,
    "cache_hits": 0
  },
  "inverted": {
    "fraction": 1.0,
    "completed_arms": 3860,
    "survivor_peaks": 3324,
    "survivor_sources": 1471,
    "survivor_dependency_groups": 48,
    "cache_hits": 0
  },
  "block_scrambled": {
    "fraction": 1.0,
    "completed_arms": 3860,
    "survivor_peaks": 803,
    "survivor_sources": 445,
    "survivor_dependency_groups": 39,
    "cache_hits": 0
  },
  "injected": {
    "fraction": 1.0,
    "completed_arms": 3860,
    "survivor_peaks": 5760,
    "survivor_sources": 2234,
    "survivor_dependency_groups": 47,
    "cache_hits": 0
  }
}
```

Injection detected/recovered: 1918/1342 of 3860 completed sources. Lead groups: 2850; holdout-supported groups: 84; unresolved: 422.

Source-level denominators are correlated; dependency groups are conservative and not certified independent. Real/control contrasts describe workload, not discovery yield. One flux-level injection per source cannot validate extraction. Host assignment, completeness, FPP and novelty remain unresolved; no promotion. Rows: reports/*.jsonl (one row per source per arm; peaks.jsonl one row per peak).
