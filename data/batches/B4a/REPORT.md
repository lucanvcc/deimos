# B4a: paired Kepler Bonus search

Status: COMPLETE. Stop: None. 5260/5260 sources terminal; 5059 complete, 201 failed, 0 interrupted (never replayed).

```json
{
  "real": {
    "fraction": 1.0,
    "completed_arms": 5059,
    "survivor_peaks": 3731,
    "survivor_sources": 1600,
    "survivor_dependency_groups": 46,
    "cache_hits": 0
  },
  "inverted": {
    "fraction": 1.0,
    "completed_arms": 5059,
    "survivor_peaks": 3853,
    "survivor_sources": 1623,
    "survivor_dependency_groups": 47,
    "cache_hits": 0
  },
  "block_scrambled": {
    "fraction": 1.0,
    "completed_arms": 5059,
    "survivor_peaks": 886,
    "survivor_sources": 509,
    "survivor_dependency_groups": 36,
    "cache_hits": 0
  },
  "injected": {
    "fraction": 1.0,
    "completed_arms": 5059,
    "survivor_peaks": 8010,
    "survivor_sources": 3039,
    "survivor_dependency_groups": 47,
    "cache_hits": 0
  }
}
```

Injection detected/recovered: 3090/2117 of 5059 completed sources. Lead groups: 3246; holdout-supported groups: 69; unresolved: 408.

Source-level denominators are correlated; dependency groups are conservative and not certified independent. Real/control contrasts describe workload, not discovery yield. One flux-level injection per source cannot validate extraction. Host assignment, completeness, FPP and novelty remain unresolved; no promotion. Rows: reports/*.jsonl (one row per source per arm; peaks.jsonl one row per peak).
