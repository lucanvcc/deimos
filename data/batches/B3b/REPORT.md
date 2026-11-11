# B3b: paired Kepler Bonus search

Status: COMPLETE. Stop: None. 4933/4933 sources terminal; 4689 complete, 244 failed, 0 interrupted (never replayed).

```json
{
  "real": {
    "fraction": 1.0,
    "completed_arms": 4689,
    "survivor_peaks": 2595,
    "survivor_sources": 1120,
    "survivor_dependency_groups": 44,
    "cache_hits": 0
  },
  "inverted": {
    "fraction": 1.0,
    "completed_arms": 4689,
    "survivor_peaks": 2518,
    "survivor_sources": 1139,
    "survivor_dependency_groups": 46,
    "cache_hits": 0
  },
  "block_scrambled": {
    "fraction": 1.0,
    "completed_arms": 4689,
    "survivor_peaks": 664,
    "survivor_sources": 399,
    "survivor_dependency_groups": 36,
    "cache_hits": 0
  },
  "injected": {
    "fraction": 1.0,
    "completed_arms": 4689,
    "survivor_peaks": 6317,
    "survivor_sources": 2563,
    "survivor_dependency_groups": 46,
    "cache_hits": 0
  }
}
```

Injection detected/recovered: 2825/1874 of 4689 completed sources. Lead groups: 2241; holdout-supported groups: 76; unresolved: 323.

Source-level denominators are correlated; dependency groups are conservative and not certified independent. Real/control contrasts describe workload, not discovery yield. One flux-level injection per source cannot validate extraction. Host assignment, completeness, FPP and novelty remain unresolved; no promotion. Rows: reports/*.jsonl (one row per source per arm; peaks.jsonl one row per peak).
