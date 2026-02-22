"""Tests for OptimalBragg.cli — command-line interface."""

import subprocess
import sys
import pytest


class TestCLIHelp:
    def test_main_help(self):
        result = subprocess.run(
            [sys.executable, '-m', 'OptimalBragg', '--help'],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        assert 'optimize' in result.stdout
        assert 'plot' in result.stdout
        assert 'mc' in result.stdout
        assert 'corner' in result.stdout

    def test_optimize_help(self):
        result = subprocess.run(
            [sys.executable, '-m', 'OptimalBragg', 'optimize', '--help'],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        assert 'params' in result.stdout

    def test_mc_help(self):
        result = subprocess.run(
            [sys.executable, '-m', 'OptimalBragg', 'mc', '--help'],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        assert 'n_samples' in result.stdout

    def test_no_command_fails(self):
        result = subprocess.run(
            [sys.executable, '-m', 'OptimalBragg'],
            capture_output=True, text=True,
        )
        assert result.returncode != 0
