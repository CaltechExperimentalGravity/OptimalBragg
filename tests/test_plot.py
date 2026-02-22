"""Tests for OptimalBragg.plot — visualization functions."""

import numpy as np
import pytest
import matplotlib
matplotlib.use('Agg')  # non-interactive backend for testing

from OptimalBragg import Material, qw_stack
from OptimalBragg.materials import SiO2, TiTa2O5, FusedSilica, air


@pytest.fixture
def aLIGO_stack():
    return qw_stack(
        lam_ref=1064e-9,
        substrate=Material(FusedSilica),
        superstrate=Material(air),
        thin_films={"L": Material(SiO2), "H": Material(TiTa2O5)},
        pattern="LH" * 8,
    )


class TestPlotStarfish:
    def test_returns_figure(self):
        from OptimalBragg.plot import plot_starfish
        costs = {'Trans1064': 3.2, 'Brownian': 1.5, 'TO': 2.1}
        fig = plot_starfish(costs)
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)

    def test_saves_to_file(self, tmp_path):
        from OptimalBragg.plot import plot_starfish
        costs = {'A': 1.0, 'B': 2.0}
        path = str(tmp_path / 'sf.png')
        fig = plot_starfish(costs, save_path=path)
        assert (tmp_path / 'sf.png').exists()
        import matplotlib.pyplot as plt
        plt.close(fig)


class TestPlotSpectral:
    def test_returns_figure(self, aLIGO_stack):
        from OptimalBragg.plot import plot_spectral
        fig = plot_spectral(
            aLIGO_stack['ns'], aLIGO_stack['Ls_opt'],
            aLIGO_stack['lam_ref'],
        )
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)


class TestPlotNoise:
    def test_returns_figure(self, aLIGO_stack):
        from OptimalBragg.plot import plot_noise
        fig = plot_noise(aLIGO_stack, w_beam=0.062)
        assert fig is not None
        import matplotlib.pyplot as plt
        plt.close(fig)


class TestSigmaClip:
    def test_removes_outliers(self):
        from OptimalBragg.plot import sigma_clip
        rng = np.random.default_rng(42)
        samples = rng.normal(0, 1, (2, 1000))
        samples[0, 0] = 100  # outlier
        clipped, mask = sigma_clip(samples, nsigma=3)
        assert clipped.shape[1] < 1000
        assert not mask[0]  # outlier removed
