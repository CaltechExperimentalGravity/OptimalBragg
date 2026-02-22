"""Micro-benchmark for multidiel1 before/after Numba JIT."""
import numpy as np
import time
from OptimalBragg.layers import multidiel1

# Build a realistic 14-bilayer stack
n = np.array([1.0] + [1.45, 3.5] * 14 + [1.45])
L = 0.25 * np.ones(28)
lamb = np.array([1.0])

# Warm up (triggers JIT compilation if Numba is active)
_ = multidiel1(n, L, lamb, 0, 'te')

# Benchmark: single wavelength
N = 10000
t0 = time.perf_counter()
for _ in range(N):
    multidiel1(n, L, lamb, 0, 'te')
dt = time.perf_counter() - t0
print(f"multidiel1 (1 wl): {1e6*dt/N:.1f} us/call  ({N} calls in {dt:.2f} s)")

# Benchmark: multi-wavelength (simulates consolidated call)
lamb3 = np.array([1.0, 0.756, 0.297])
t0 = time.perf_counter()
for _ in range(N):
    multidiel1(n, L, lamb3, 0, 'te')
dt = time.perf_counter() - t0
print(f"multidiel1 (3 wl): {1e6*dt/N:.1f} us/call  ({N} calls in {dt:.2f} s)")
