# B3f: paired Kepler Bonus search

Status: COMPLETE. Stop: None. 5000/5000 sources terminal; 4614 complete, 386 failed, 0 interrupted (never replayed).

```json
{
  "real": {
    "fraction": 1.0,
    "completed_arms": 4614,
    "survivor_peaks": 2553,
    "survivor_sources": 1119,
    "survivor_dependency_groups": 47,
    "cache_hits": 0
  },
  "inverted": {
    "fraction": 1.0,
    "completed_arms": 4614,
    "survivor_peaks": 2820,
    "survivor_sources": 1270,
    "survivor_dependency_groups": 48,
    "cache_hits": 0
  },
  "block_scrambled": {
    "fraction": 1.0,
    "completed_arms": 4614,
    "survivor_peaks": 615,
    "survivor_sources": 351,
    "survivor_dependency_groups": 40,
    "cache_hits": 0
  },
  "injected": {
    "fraction": 1.0,
    "completed_arms": 4614,
    "survivor_peaks": 6620,
    "survivor_sources": 2594,
    "survivor_dependency_groups": 48,
    "cache_hits": 0
  }
}
```

Injection detected/recovered: 2728/1925 of 4614 completed sources. Lead groups: 2342; holdout-supported groups: 77; unresolved: 366.

Source-level denominators are correlated; dependency groups are conservative and not certified independent. Real/control contrasts describe workload, not discovery yield. One flux-level injection per source cannot validate extraction. Host assignment, completeness, FPP and novelty remain unresolved; no promotion. Rows: reports/*.jsonl (one row per source per arm; peaks.jsonl one row per peak).
