"""Benchmark: JIT vs pure-numpy thermooptic noise calculation."""
import time
import numpy as np

from OptimalBragg import Material, qw_stack
from OptimalBragg.materials import SiO2, TiTa2O5, FusedSilica, air
from OptimalBragg.noise import (
    coating_thermooptic_fast, coating_thermooptic,
    extract_stack_params,
)

stack = qw_stack(
    lam_ref=1064e-9,
    substrate=Material(FusedSilica),
    superstrate=Material(air),
    thin_films={"L": Material(SiO2), "H": Material(TiTa2O5)},
    pattern="LH" * 14,
)
stack_params = extract_stack_params(stack, 0.17, 0.20)

L = stack["Ls_opt"]
f = 100.0
w_beam = 0.062
N = 10_000


def bench_numpy():
    # warm up
    coating_thermooptic(np.array([f]), stack, w_beam, 0.17, 0.20)

    t0 = time.perf_counter()
    for _ in range(N):
        coating_thermooptic(np.array([f]), stack, w_beam, 0.17, 0.20)
    dt = time.perf_counter() - t0
    return dt


def bench_jit():
    # warm up JIT
    coating_thermooptic_fast(f, L, 1064e-9, w_beam, stack_params)

    t0 = time.perf_counter()
    for _ in range(N):
        coating_thermooptic_fast(f, L, 1064e-9, w_beam, stack_params)
    dt = time.perf_counter() - t0
    return dt


if __name__ == '__main__':
    n_layers = len(L)
    print(f"Benchmarking thermooptic noise ({N} calls, {n_layers}-layer "
          f"stack, f=100 Hz)")
    print()

    dt_np = bench_numpy()
    us_np = 1e6 * dt_np / N
    print(f"  numpy:  {us_np:.1f} us/call  ({N} calls in {dt_np:.2f}s)")

    dt_jit = bench_jit()
    us_jit = 1e6 * dt_jit / N
    print(f"  JIT:    {us_jit:.1f} us/call  ({N} calls in {dt_jit:.2f}s)")

    print(f"\n  Speedup: {us_np/us_jit:.1f}x")

    # Verify agreement
    StoZ_np, _, _ = coating_thermooptic(
        np.array([f]), stack, w_beam, 0.17, 0.20)
    StoZ_jit = coating_thermooptic_fast(
        f, L, 1064e-9, w_beam, stack_params)
    print(f"\n  numpy result: {StoZ_np:.15e}")
    print(f"  JIT result:   {StoZ_jit:.15e}")
    print(f"  Agreement within 5%: {np.isclose(StoZ_jit, StoZ_np, rtol=0.05)}")
