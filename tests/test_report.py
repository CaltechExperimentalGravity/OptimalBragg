"""Tests for OptimalBragg.report — Sphinx report generation."""

import pytest


class TestAssessQuality:
    def test_pass_close(self):
        from OptimalBragg.report import assess_quality
        assert assess_quality(5e-6, 5e-6, 'Trans1') == 'PASS'

    def test_warn_moderate(self):
        from OptimalBragg.report import assess_quality
        assert assess_quality(8e-6, 5e-6, 'Trans1') == 'WARN'

    def test_fail_far(self):
        from OptimalBragg.report import assess_quality
        assert assess_quality(50e-6, 5e-6, 'Trans1') == 'FAIL'

    def test_zero_target(self):
        from OptimalBragg.report import assess_quality
        assert assess_quality(1.0, 0, 'Brownian') == 'PASS'

    def test_brownian_loose(self):
        from OptimalBragg.report import assess_quality
        # 50% off is OK for non-transmission costs
        assert assess_quality(15, 20, 'Brownian') == 'PASS'


class TestFormatValue:
    def test_trans1_hr(self):
        """HR target (T=5 ppm) should format as ppm."""
        from OptimalBragg.report import _format_value
        result = _format_value(5e-6, 'Trans1', target=5e-6)
        assert 'ppm' in result

    def test_trans2_ar(self):
        """AR target (T=1.0) should format as reflectivity ppm."""
        from OptimalBragg.report import _format_value
        result = _format_value(0.95, 'Trans2', target=1.0)
        assert 'ppm' in result

    def test_trans2_hr(self):
        """HR target (T=0.032) should format as ppm."""
        from OptimalBragg.report import _format_value
        result = _format_value(0.032, 'Trans2', target=0.032)
        assert 'ppm' in result


class TestTransLabel:
    def test_hr_label(self):
        from OptimalBragg.report import _trans_label
        label = _trans_label('Trans1', 5e-6, 1064)
        assert label == 'T @ 1064 nm'

    def test_ar_label(self):
        from OptimalBragg.report import _trans_label
        label = _trans_label('Trans2', 1.0, 700)
        assert label == 'R @ 700 nm'
