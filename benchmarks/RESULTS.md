# Benchmark Results

Run `python benchmarks/bench_multidiel1.py` from the repo root to reproduce.

## multidiel1 micro-benchmark (14-bilayer stack, 10k calls)

| Date | Stage | multidiel1 (1 wl) | multidiel1 (3 wl) | Notes |
|------|-------|-------------------|-------------------|-------|
| 2026-02-18 | After Numba JIT + call consolidation | 7.4 us/call | 7.8 us/call | 3 wl barely slower than 1 wl |

## Full optimization (SiN_aSi mkETM, 14 pairs, popsize=200)

| Date | Wall time | Iterations | Final cost | TransPSL | Brownian | Thermooptic |
|------|-----------|------------|------------|----------|----------|-------------|
| 2026-02-18 | 31.0 sec | 87 + polish | 1.919 | 0.0001 | 0.6721 | 0.2838 |
