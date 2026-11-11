# B3e: paired Kepler Bonus search

Status: COMPLETE. Stop: None. 5000/5000 sources terminal; 4632 complete, 368 failed, 0 interrupted (never replayed).

```json
{
  "real": {
    "fraction": 1.0,
    "completed_arms": 4632,
    "survivor_peaks": 2931,
    "survivor_sources": 1270,
    "survivor_dependency_groups": 47,
    "cache_hits": 0
  },
  "inverted": {
    "fraction": 1.0,
    "completed_arms": 4632,
    "survivor_peaks": 3096,
    "survivor_sources": 1315,
    "survivor_dependency_groups": 46,
    "cache_hits": 0
  },
  "block_scrambled": {
    "fraction": 1.0,
    "completed_arms": 4632,
    "survivor_peaks": 830,
    "survivor_sources": 448,
    "survivor_dependency_groups": 41,
    "cache_hits": 0
  },
  "injected": {
    "fraction": 1.0,
    "completed_arms": 4632,
    "survivor_peaks": 6668,
    "survivor_sources": 2600,
    "survivor_dependency_groups": 49,
    "cache_hits": 0
  }
}
```

Injection detected/recovered: 2654/1815 of 4632 completed sources. Lead groups: 2638; holdout-supported groups: 86; unresolved: 404.

Source-level denominators are correlated; dependency groups are conservative and not certified independent. Real/control contrasts describe workload, not discovery yield. One flux-level injection per source cannot validate extraction. Host assignment, completeness, FPP and novelty remain unresolved; no promotion. Rows: reports/*.jsonl (one row per source per arm; peaks.jsonl one row per peak).
