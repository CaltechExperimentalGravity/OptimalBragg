"""Tests for OptimalBragg.report — Sphinx report generation."""

import pytest


class TestAssessQuality:
    def test_pass_close(self):
        from OptimalBragg.report import assess_quality
        assert assess_quality(5e-6, 5e-6, 'Trans1064') == 'PASS'

    def test_warn_moderate(self):
        from OptimalBragg.report import assess_quality
        assert assess_quality(8e-6, 5e-6, 'Trans1064') == 'WARN'

    def test_fail_far(self):
        from OptimalBragg.report import assess_quality
        assert assess_quality(50e-6, 5e-6, 'Trans1064') == 'FAIL'

    def test_zero_target(self):
        from OptimalBragg.report import assess_quality
        assert assess_quality(1.0, 0, 'Brownian') == 'PASS'

    def test_brownian_loose(self):
        from OptimalBragg.report import assess_quality
        # 50% off is OK for non-transmission costs
        assert assess_quality(15, 20, 'Brownian') == 'PASS'


class TestFormatValue:
    def test_trans1064(self):
        from OptimalBragg.report import _format_value
        result = _format_value(5e-6, 'Trans1064')
        assert 'ppm' in result

    def test_trans532(self):
        from OptimalBragg.report import _format_value
        result = _format_value(0.032, 'Trans532')
        assert '%' in result
