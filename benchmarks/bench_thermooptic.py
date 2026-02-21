"""Benchmark: gwinc vs JIT thermooptic noise calculation."""
import copy
import time
import numpy as np
import gwinc
import gwinc.noise.coatingthermal

from generic.thermoopticUtils import coating_thermooptic_fast, extract_ifo_params

ifo = gwinc.Struct.from_file('SiN_aSi/aSiSiN.yaml')
ifo_params = extract_ifo_params(ifo)

L = 0.25 * np.ones(29)
f = 100.0
N = 10_000


def bench_gwinc():
    mir = copy.copy(ifo.Optics.ETM)
    mir.Coating = copy.copy(ifo.Optics.ETM.Coating)
    mir.Coating.dOpt = L
    # warm up
    gwinc.noise.coatingthermal.coating_thermooptic(
        f, mir, ifo.Laser.Wavelength, ifo.Optics.ETM.BeamRadius)

    t0 = time.perf_counter()
    for _ in range(N):
        mir.Coating.dOpt = L
        gwinc.noise.coatingthermal.coating_thermooptic(
            f, mir, ifo.Laser.Wavelength, ifo.Optics.ETM.BeamRadius)
    dt = time.perf_counter() - t0
    return dt


def bench_jit():
    # warm up JIT
    coating_thermooptic_fast(f, L, ifo.Laser.Wavelength,
                             ifo.Optics.ETM.BeamRadius, ifo_params)

    t0 = time.perf_counter()
    for _ in range(N):
        coating_thermooptic_fast(f, L, ifo.Laser.Wavelength,
                                 ifo.Optics.ETM.BeamRadius, ifo_params)
    dt = time.perf_counter() - t0
    return dt


if __name__ == '__main__':
    print(f"Benchmarking thermooptic noise ({N} calls, 29-layer stack, f=100 Hz)")
    print()

    dt_gwinc = bench_gwinc()
    us_gwinc = 1e6 * dt_gwinc / N
    print(f"  gwinc:  {us_gwinc:.1f} us/call  ({N} calls in {dt_gwinc:.2f}s)")

    dt_jit = bench_jit()
    us_jit = 1e6 * dt_jit / N
    print(f"  JIT:    {us_jit:.1f} us/call  ({N} calls in {dt_jit:.2f}s)")

    print(f"\n  Speedup: {us_gwinc/us_jit:.1f}x")

    # Verify correctness
    mir = copy.copy(ifo.Optics.ETM)
    mir.Coating = copy.copy(ifo.Optics.ETM.Coating)
    mir.Coating.dOpt = L
    StoZ_gwinc = gwinc.noise.coatingthermal.coating_thermooptic(
        f, mir, ifo.Laser.Wavelength, ifo.Optics.ETM.BeamRadius)[0]
    StoZ_jit = coating_thermooptic_fast(
        f, L, ifo.Laser.Wavelength, ifo.Optics.ETM.BeamRadius, ifo_params)
    print(f"\n  gwinc result: {StoZ_gwinc:.15e}")
    print(f"  JIT result:   {StoZ_jit:.15e}")
    print(f"  Match: {np.isclose(StoZ_jit, StoZ_gwinc, rtol=1e-12)}")
