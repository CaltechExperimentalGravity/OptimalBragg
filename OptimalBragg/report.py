"""Sphinx run report generation utilities.

Generates .rst files for individual optimization runs with inline plots,
summary tables, and quality assessment. Also maintains an auto-generated
index of all run reports.
"""

import os
import re
import subprocess
import h5py
import numpy as np
from datetime import datetime
from pathlib import Path

# Project root: parent of OptimalBragg/
_PROJECT_ROOT = str(Path(__file__).parent.parent)


def assess_quality(achieved, target, metric_name):
    """Compare an achieved value to its target.

    Returns
    -------
    status : str
        'PASS', 'WARN', or 'FAIL'.
    """
    if target == 0:
        return 'PASS'
    ratio = abs(achieved - target) / abs(target)
    if metric_name.startswith('Trans'):
        if ratio < 0.5:
            return 'PASS'
        elif ratio < 2.0:
            return 'WARN'
        else:
            return 'FAIL'
    if ratio < 1.0:
        return 'PASS'
    elif ratio < 5.0:
        return 'WARN'
    return 'FAIL'


def _format_value(val, metric_name):
    """Format a metric value for display in the summary table."""
    if metric_name in ('Trans1064',):
        return f'{val*1e6:.2f} ppm'
    elif metric_name in ('Trans532', 'TransOPLEV'):
        return f'{val*100:.3f}%'
    elif metric_name == 'Brownian':
        return f'{val:.4f}'
    elif metric_name == 'Thermooptic':
        return f'{val:.4e}'
    return f'{val:.4g}'


def _format_target(target, metric_name):
    """Format a target value for display."""
    if metric_name in ('Trans1064',):
        return f'{target*1e6:.1f} ppm'
    elif metric_name in ('Trans532', 'TransOPLEV'):
        return f'{target*100:.1f}%'
    elif metric_name == 'Brownian':
        return f'{target:.2f}'
    elif metric_name == 'Thermooptic':
        return f'{target:.1e}'
    return f'{target:.4g}'


def generate_run_rst(hdf5_path, mirror_type, fig_paths, params,
                     project_dir='Arms', mc_hdf5_path=None):
    """Generate a .rst file for one optimization run.

    Parameters
    ----------
    hdf5_path : str
        Path to the optimizer output HDF5 file.
    mirror_type : str
        'ETM' or 'ITM'.
    fig_paths : dict
        Mapping of plot names to file paths.
    params : dict
        Optimization parameters (from YAML).
    project_dir : str
        Project directory name.
    mc_hdf5_path : str, optional
        Path to Monte Carlo output HDF5.

    Returns
    -------
    rst_content : str
    rst_path : str
    """
    from OptimalBragg.layers import op2phys

    # Extract timestamp from HDF5 filename
    basename = os.path.basename(hdf5_path)
    m = re.search(r'(\d{6}_\d{6})', basename)
    timestamp = m.group(1) if m else datetime.now().strftime('%y%m%d_%H%M%S')

    # Read metrics from HDF5
    metrics = {}
    with h5py.File(hdf5_path, 'r') as f:
        grp = 'diffevo_output'
        for key in ('T1064', 'T532', 'TOPL', 'scalarCost'):
            try:
                metrics[key] = float(np.array(f[f'{grp}/{key}']))
            except (KeyError, TypeError):
                pass
        try:
            L = np.array(f[f'{grp}/L'])
            n = np.array(f[f'{grp}/n'])
            metrics['nLayers'] = len(L)
        except KeyError:
            L, n = None, None

        vector_costs = {}
        try:
            for cost_name in f[f'{grp}/vectorCost']:
                vector_costs[cost_name] = float(
                    np.array(f[f'{grp}/vectorCost/{cost_name}']))
        except KeyError:
            pass

    # Compute physical thickness
    if L is not None and n is not None:
        try:
            wl = params.get('_wavelength', 1064e-9)
            L_phys = wl * op2phys(L, n[1:-1])
            metrics['thickness_um'] = np.sum(L_phys) * 1e6
        except Exception:
            pass

    # Build .rst
    title = f'{mirror_type} Run {timestamp}'
    lines = [title, '=' * len(title), '']

    # Summary table
    lines.extend([
        '.. list-table:: Design Summary',
        '   :header-rows: 1',
        '   :widths: 30 25 25 20',
        '',
        '   * - Metric',
        '     - Achieved',
        '     - Target',
        '     - Status',
    ])

    metric_to_cost = {
        'T1064': 'Trans1064',
        'T532': 'Trans532',
        'TOPL': 'TransOPLEV',
    }

    for hdf_key, cost_name in metric_to_cost.items():
        if hdf_key in metrics and cost_name in params.get('costs', {}):
            achieved = metrics[hdf_key]
            target = params['costs'][cost_name]['target']
            weight = params['costs'][cost_name]['weight']
            if weight == 0:
                continue
            status = assess_quality(achieved, target, cost_name)
            lines.extend([
                f'   * - {cost_name}',
                f'     - {_format_value(achieved, cost_name)}',
                f'     - {_format_target(target, cost_name)}',
                f'     - {status}',
            ])

    if 'nLayers' in metrics:
        lines.extend([
            f'   * - Number of layers',
            f'     - {metrics["nLayers"]}',
            f'     -',
            f'     -',
        ])
    if 'thickness_um' in metrics:
        lines.extend([
            f'   * - Total thickness',
            f'     - {metrics["thickness_um"]:.2f} um',
            f'     -',
            f'     -',
        ])
    if 'scalarCost' in metrics:
        lines.extend([
            f'   * - Total scalar cost',
            f'     - {metrics["scalarCost"]:.6f}',
            f'     -',
            f'     -',
        ])
    lines.append('')

    # Cost breakdown
    if vector_costs:
        lines.extend([
            '.. list-table:: Cost Function Breakdown',
            '   :header-rows: 1',
            '   :widths: 40 30 30',
            '',
            '   * - Cost Term',
            '     - Value',
            '     - Weight',
        ])
        for cost_name, cost_val in vector_costs.items():
            weight = params.get('costs', {}).get(
                cost_name, {}).get('weight', 0)
            lines.extend([
                f'   * - {cost_name}',
                f'     - {cost_val:.6f}',
                f'     - {weight}',
            ])
        lines.append('')

    # Plot sections
    rel_prefix = '../../..'

    if 'layers' in fig_paths:
        lines.extend([
            'Layer Structure & E-Field',
            '-------------------------', '',
            f'.. image:: {rel_prefix}/{fig_paths["layers"]}',
            '   :width: 100%', '',
        ])

    if 'spectral' in fig_paths:
        lines.extend([
            'Spectral Reflectivity',
            '---------------------', '',
            f'.. image:: {rel_prefix}/{fig_paths["spectral"]}',
            '   :width: 100%', '',
        ])

    if 'starfish' in fig_paths:
        lines.extend([
            'Starfish Cost Chart',
            '-------------------', '',
            f'.. image:: {rel_prefix}/{fig_paths["starfish"]}',
            '   :width: 60%', '',
        ])

    if 'thermal_noise' in fig_paths:
        lines.extend([
            'Thermal Noise Budget',
            '--------------------', '',
            f'.. image:: {rel_prefix}/{fig_paths["thermal_noise"]}',
            '   :width: 100%', '',
        ])

    # MC section
    mc_labels = [
        ('T1064_ppm', r':math:`T_{1064}` [ppm]'),
        ('T532_pct', r':math:`T_{532}` [%]'),
        ('S_TO', r':math:`S_\mathrm{TO}` '
         r'[:math:`\times 10^{-21}` m/:math:`\sqrt{\mathrm{Hz}}`]'),
        ('S_Br', r':math:`S_\mathrm{Br}` '
         r'[:math:`\times 10^{-21}` m/:math:`\sqrt{\mathrm{Hz}}`]'),
        ('E_surf', r':math:`E_\mathrm{surface}` [V/m]'),
    ]

    if mc_hdf5_path and os.path.isfile(mc_hdf5_path):
        lines.extend([
            'Monte Carlo Sensitivity',
            '-----------------------', '',
        ])
        try:
            with h5py.File(mc_hdf5_path, 'r') as fmc:
                mc_samples = np.array(fmc['MCout'][:])
            n_obs, n_mc = mc_samples.shape
            lines.extend([
                f'Ensemble of {n_mc} MC samples with 0.5% Gaussian '
                f'perturbations on refractive indices and layer thicknesses.',
                '',
                '.. list-table:: MC Statistics (median, 5th-95th percentile)',
                '   :header-rows: 1',
                '   :widths: 40 20 20 20', '',
                '   * - Observable',
                '     - Median',
                '     - 5th %ile',
                '     - 95th %ile',
            ])
            for i in range(min(n_obs, len(mc_labels))):
                row = mc_samples[i]
                med = np.median(row)
                lo = np.percentile(row, 5)
                hi = np.percentile(row, 95)
                _, label = mc_labels[i]
                lines.extend([
                    f'   * - {label}',
                    f'     - {med:.4g}',
                    f'     - {lo:.4g}',
                    f'     - {hi:.4g}',
                ])
            lines.append('')
        except Exception:
            lines.extend(['*Could not read MC data.*', ''])

        if 'corner' in fig_paths:
            lines.extend([
                f'.. image:: {rel_prefix}/{fig_paths["corner"]}',
                '   :width: 100%', '',
            ])
    else:
        lines.extend([
            'Monte Carlo Sensitivity',
            '-----------------------', '',
            '*MC analysis pending. Re-run report generation after '
            'MC completes.*', '',
        ])

    rst_content = '\n'.join(lines) + '\n'

    # Write .rst
    rst_dir = os.path.join(_PROJECT_ROOT, 'docs', 'runs', mirror_type)
    os.makedirs(rst_dir, exist_ok=True)
    rst_path = os.path.join(rst_dir, f'{timestamp}.rst')
    with open(rst_path, 'w') as f:
        f.write(rst_content)

    # Regenerate index and rebuild HTML
    regenerate_runs_index(os.path.join(_PROJECT_ROOT, 'docs', 'runs'))
    build_html()

    return rst_content, rst_path


def regenerate_runs_index(runs_dir):
    """Scan docs/runs/ for all .rst files, rebuild index.rst."""
    os.makedirs(runs_dir, exist_ok=True)

    entries = {}
    for mirror_type in ('ETM', 'ITM'):
        subdir = os.path.join(runs_dir, mirror_type)
        if not os.path.isdir(subdir):
            continue
        timestamps = []
        for fname in sorted(os.listdir(subdir), reverse=True):
            if fname.endswith('.rst') and fname != 'index.rst':
                timestamps.append(fname[:-4])
        if timestamps:
            entries[mirror_type] = timestamps

    lines = [
        'Run Reports', '===========', '',
        'Auto-generated index of optimization run reports.', '',
    ]

    for mirror_type in ('ETM', 'ITM'):
        if mirror_type not in entries:
            continue
        lines.extend([
            f'{mirror_type} Runs',
            '-' * (len(mirror_type) + 5), '',
            '.. toctree::', '   :maxdepth: 1', '',
        ])
        for ts in entries[mirror_type]:
            lines.append(f'   {mirror_type}/{ts}')
        lines.append('')

    if not entries:
        lines.extend([
            'No run reports found yet.', '',
        ])

    index_path = os.path.join(runs_dir, 'index.rst')
    with open(index_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def build_html():
    """Run ``make html`` in docs/ to rebuild Sphinx pages."""
    docs_dir = os.path.join(_PROJECT_ROOT, 'docs')
    if not os.path.isfile(os.path.join(docs_dir, 'Makefile')):
        return
    try:
        result = subprocess.run(
            ['make', 'html'], cwd=docs_dir,
            capture_output=True, text=True, timeout=120,
        )
        if result.returncode == 0:
            print(f"Sphinx HTML built: {docs_dir}/_build/html/")
        else:
            print(f"Sphinx build warning (exit {result.returncode}): "
                  f"{result.stderr[-200:] if result.stderr else ''}")
    except FileNotFoundError:
        print("Sphinx build skipped: 'make' not found")
    except subprocess.TimeoutExpired:
        print("Sphinx build skipped: timed out after 120s")
