# B3a: paired Kepler Bonus search

Status: COMPLETE. Stop: None. 1000/1000 sources terminal; 980 complete, 20 failed, 0 interrupted (never replayed).

```json
{
  "real": {
    "fraction": 1.0,
    "completed_arms": 980,
    "survivor_peaks": 572,
    "survivor_sources": 249,
    "survivor_dependency_groups": 32,
    "cache_hits": 0
  },
  "inverted": {
    "fraction": 1.0,
    "completed_arms": 980,
    "survivor_peaks": 592,
    "survivor_sources": 274,
    "survivor_dependency_groups": 35,
    "cache_hits": 0
  },
  "block_scrambled": {
    "fraction": 1.0,
    "completed_arms": 980,
    "survivor_peaks": 183,
    "survivor_sources": 105,
    "survivor_dependency_groups": 22,
    "cache_hits": 0
  },
  "injected": {
    "fraction": 1.0,
    "completed_arms": 980,
    "survivor_peaks": 1337,
    "survivor_sources": 537,
    "survivor_dependency_groups": 34,
    "cache_hits": 0
  }
}
```

Injection detected/recovered: 592/385 of 980 completed sources. Lead groups: 542; holdout-supported groups: 8; unresolved: 38.

Source-level denominators are correlated; dependency groups are conservative and not certified independent. Real/control contrasts describe workload, not discovery yield. One flux-level injection per source cannot validate extraction. Host assignment, completeness, FPP and novelty remain unresolved; no promotion. Rows: reports/*.jsonl (one row per source per arm; peaks.jsonl one row per peak).
