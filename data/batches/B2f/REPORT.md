# B2f: paired Kepler Bonus search

Status: COMPLETE. Stop: None. 1634/1634 sources terminal; 1572 complete, 55 failed, 7 interrupted (never replayed).

```json
{
  "real": {
    "fraction": 1.0,
    "completed_arms": 1572,
    "survivor_peaks": 1645,
    "survivor_sources": 687,
    "survivor_dependency_groups": 42,
    "cache_hits": 0
  },
  "inverted": {
    "fraction": 1.0,
    "completed_arms": 1572,
    "survivor_peaks": 1715,
    "survivor_sources": 731,
    "survivor_dependency_groups": 43,
    "cache_hits": 0
  },
  "block_scrambled": {
    "fraction": 1.0,
    "completed_arms": 1572,
    "survivor_peaks": 255,
    "survivor_sources": 146,
    "survivor_dependency_groups": 35,
    "cache_hits": 0
  },
  "injected": {
    "fraction": 1.0,
    "completed_arms": 1572,
    "survivor_peaks": 3897,
    "survivor_sources": 1286,
    "survivor_dependency_groups": 43,
    "cache_hits": 0
  }
}
```

Injection detected/recovered: 1262/980 of 1572 completed sources. Lead groups: 1488; holdout-supported groups: 22; unresolved: 186.

Source-level denominators are correlated; dependency groups are conservative and not certified independent. Real/control contrasts describe workload, not discovery yield. One flux-level injection per source cannot validate extraction. Host assignment, completeness, FPP and novelty remain unresolved; no promotion. Rows: reports/*.jsonl (one row per source per arm; peaks.jsonl one row per peak).
