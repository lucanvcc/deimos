# B3k: paired Kepler Bonus search

Status: COMPLETE. Stop: None. 145/145 sources terminal; 86 complete, 59 failed, 0 interrupted (never replayed).

```json
{
  "real": {
    "fraction": 1.0,
    "completed_arms": 86,
    "survivor_peaks": 86,
    "survivor_sources": 38,
    "survivor_dependency_groups": 15,
    "cache_hits": 0
  },
  "inverted": {
    "fraction": 1.0,
    "completed_arms": 86,
    "survivor_peaks": 123,
    "survivor_sources": 50,
    "survivor_dependency_groups": 13,
    "cache_hits": 0
  },
  "block_scrambled": {
    "fraction": 1.0,
    "completed_arms": 86,
    "survivor_peaks": 12,
    "survivor_sources": 4,
    "survivor_dependency_groups": 4,
    "cache_hits": 0
  },
  "injected": {
    "fraction": 1.0,
    "completed_arms": 86,
    "survivor_peaks": 178,
    "survivor_sources": 64,
    "survivor_dependency_groups": 17,
    "cache_hits": 0
  }
}
```

Injection detected/recovered: 57/45 of 86 completed sources. Lead groups: 83; holdout-supported groups: 7; unresolved: 7.

Source-level denominators are correlated; dependency groups are conservative and not certified independent. Real/control contrasts describe workload, not discovery yield. One flux-level injection per source cannot validate extraction. Host assignment, completeness, FPP and novelty remain unresolved; no promotion. Rows: reports/*.jsonl (one row per source per arm; peaks.jsonl one row per peak).
