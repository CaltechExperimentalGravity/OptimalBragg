"""Optimize ETM coating for LIGO Voyager (aSi/SiN at 2050 nm, 123 K).

Usage:
    cd projects/Voyager_aSiSiN
    python mkMirror.py
"""
from OptimalBragg.optimizer import run_optimization

if __name__ == '__main__':
    run_optimization('ETM_params.yml')
