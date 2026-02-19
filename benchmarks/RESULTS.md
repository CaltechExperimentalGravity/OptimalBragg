# Benchmark Results

Run `python benchmarks/bench_multidiel1.py` from the repo root to reproduce.

| Date | Stage | multidiel1 (1 wl) | multidiel1 (3 wl) | Notes |
|------|-------|-------------------|-------------------|-------|
| *pending* | Baseline (pre-refactor) | — | — | Run before any code changes |
| *pending* | After Numba JIT | — | — | Expected 10-50x speedup |
| *pending* | After call consolidation | — | — | 3x fewer calls per getMirrorCost |
