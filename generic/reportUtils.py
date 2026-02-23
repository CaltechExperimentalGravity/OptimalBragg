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

# Project root: parent of generic/
_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def assess_quality(achieved, target, metric_name):
    """Compare an achieved value to its target.

    Parameters
    ----------
    achieved : float
        Achieved value from the optimizer.
    target : float
        Target value from the YAML config.
    metric_name : str
        Name of the metric (used for tolerance selection).

    Returns
    -------
    status : str
        'PASS', 'WARN', or 'FAIL'.
    """
    if target == 0:
        return 'PASS'
    ratio = abs(achieved - target) / abs(target)
    # Transmission targets: tight tolerance
    if metric_name.startswith('Trans'):
        if ratio < 0.5:
            return 'PASS'
        elif ratio < 2.0:
            return 'WARN'
        else:
            return 'FAIL'
    # Other costs: looser tolerance
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
        Mapping of plot names to file paths, e.g.
        ``{'layers': 'Figures/ETM/ETM_Layers_260220_073600.pdf', ...}``.
    params : dict
        Optimization parameters (from importParams).
    project_dir : str
        Project directory name ('Arms' or 'SiN_aSi').
    mc_hdf5_path : str, optional
        Path to Monte Carlo output HDF5 (with ``MCout`` dataset).
        If provided, MC statistics and corner plot are included.

    Returns
    -------
    rst_content : str
        Complete .rst page content.
    rst_path : str
        Path where the .rst file was written.
    """
    # Extract timestamp from HDF5 filename
    basename = os.path.basename(hdf5_path)
    # Match pattern like ETM_Layers_260220_073600.hdf5
    m = re.search(r'(\d{6}_\d{6})', basename)
    if m:
        timestamp = m.group(1)
    else:
        timestamp = datetime.now().strftime('%y%m%d_%H%M%S')

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

        # Vector costs
        vector_costs = {}
        try:
            for cost_name in f[f'{grp}/vectorCost']:
                vector_costs[cost_name] = float(
                    np.array(f[f'{grp}/vectorCost/{cost_name}']))
        except KeyError:
            pass

    # Compute physical thickness if possible
    if L is not None and n is not None:
        from generic.coatingUtils import op2phys
        try:
            wl = params.get('_wavelength', 1064e-9)
            L_phys = wl * op2phys(L, n[1:-1])
            metrics['thickness_um'] = np.sum(L_phys) * 1e6
        except Exception:
            pass

    # Build the .rst content
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

    # Map HDF5 keys to cost names
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

    # Individual cost term values
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
            weight = params.get('costs', {}).get(cost_name, {}).get('weight', 0)
            lines.extend([
                f'   * - {cost_name}',
                f'     - {cost_val:.6f}',
                f'     - {weight}',
            ])
        lines.append('')

    # Plot sections
    # Relative paths from docs/runs/{ETM,ITM}/ to project figures
    # From docs/runs/{ETM,ITM}/ up 3 levels to project root
    rel_prefix = f'../../../{project_dir}'

    if 'layers' in fig_paths:
        lines.extend([
            'Layer Structure & E-Field',
            '-------------------------',
            '',
            f'.. image:: {rel_prefix}/{fig_paths["layers"]}',
            '   :width: 100%',
            '',
            'Electric field profile through the dielectric stack. The upper panel',
            r'shows the normalized :math:`|E(z)|^2`, and the lower panel shows the physical',
            'layer thicknesses. Surface E-field and integrated absorption are annotated.',
            '',
        ])

    if 'reflectivity' in fig_paths:
        lines.extend([
            'Spectral Reflectivity',
            '---------------------',
            '',
            f'.. image:: {rel_prefix}/{fig_paths["reflectivity"]}',
            '   :width: 100%',
            '',
            'Spectral transmission and reflectivity across the operating band.',
            'The inset starfish chart shows the relative quality of each cost',
            'function term (higher = better).',
            '',
        ])

    if 'starfish' in fig_paths:
        lines.extend([
            'Starfish Cost Chart',
            '-------------------',
            '',
            f'.. image:: {rel_prefix}/{fig_paths["starfish"]}',
            '   :width: 60%',
            '',
            'Polar chart of cost function term quality. Larger radii indicate',
            'better performance on that objective.',
            '',
        ])

    if 'thermal_noise' in fig_paths:
        lines.extend([
            'Thermal Noise Budget',
            '--------------------',
            '',
            f'.. image:: {rel_prefix}/{fig_paths["thermal_noise"]}',
            '   :width: 100%',
            '',
            'Displacement noise budget showing Brownian, thermo-optic',
            '(thermo-elastic + thermo-refractive), and substrate contributions.',
            '',
        ])

    # Monte Carlo sensitivity section
    mc_labels = [
        ('T1064_ppm', r':math:`T_{1064}` [ppm]'),
        ('T532_pct', r':math:`T_{532}` [%]'),
        ('S_TO', r':math:`S_\mathrm{TO}` [:math:`\times 10^{-21}` m/:math:`\sqrt{\mathrm{Hz}}`]'),
        ('S_Br', r':math:`S_\mathrm{Br}` [:math:`\times 10^{-21}` m/:math:`\sqrt{\mathrm{Hz}}`]'),
        ('E_surf', r':math:`E_\mathrm{surface}` [V/m]'),
    ]

    if mc_hdf5_path and os.path.isfile(mc_hdf5_path):
        lines.extend([
            'Monte Carlo Sensitivity',
            '-----------------------',
            '',
        ])
        try:
            with h5py.File(mc_hdf5_path, 'r') as fmc:
                mc_samples = np.array(fmc['MCout'][:])
            n_obs, n_mc = mc_samples.shape
            lines.extend([
                f'Ensemble of {n_mc} Monte Carlo samples with 0.5% Gaussian perturbations',
                f'on angle of incidence, refractive indices, and layer thicknesses.',
                '',
                '.. list-table:: MC Statistics (median, 5th-95th percentile)',
                '   :header-rows: 1',
                '   :widths: 40 20 20 20',
                '',
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
                '   :width: 100%',
                '',
                'Corner plot showing pairwise correlations and marginal distributions',
                'from the Monte Carlo sensitivity analysis.',
                '',
            ])
    else:
        lines.extend([
            'Monte Carlo Sensitivity',
            '-----------------------',
            '',
            '*Monte Carlo analysis pending. Re-run report generation after MC completes.*',
            '',
        ])

    rst_content = '\n'.join(lines) + '\n'

    # Write the .rst file (always relative to project root)
    rst_dir = os.path.join(_PROJECT_ROOT, 'docs', 'runs', mirror_type)
    os.makedirs(rst_dir, exist_ok=True)
    rst_path = os.path.join(rst_dir, f'{timestamp}.rst')
    with open(rst_path, 'w') as f:
        f.write(rst_content)

    # Regenerate the index and rebuild HTML
    regenerate_runs_index(os.path.join(_PROJECT_ROOT, 'docs', 'runs'))
    build_html()

    return rst_content, rst_path


def regenerate_runs_index(runs_dir):
    """Scan docs/runs/ for all .rst files, rebuild index.rst.

    Parameters
    ----------
    runs_dir : str
        Path to the docs/runs/ directory.
    """
    os.makedirs(runs_dir, exist_ok=True)

    entries = {}  # {mirror_type: [timestamp, ...]}
    for mirror_type in ('ETM', 'ITM'):
        subdir = os.path.join(runs_dir, mirror_type)
        if not os.path.isdir(subdir):
            continue
        timestamps = []
        for fname in sorted(os.listdir(subdir), reverse=True):
            if fname.endswith('.rst') and fname != 'index.rst':
                timestamps.append(fname[:-4])  # strip .rst
        if timestamps:
            entries[mirror_type] = timestamps

    lines = [
        'Run Reports',
        '===========',
        '',
        'Auto-generated index of optimization run reports.',
        '',
    ]

    for mirror_type in ('ETM', 'ITM'):
        if mirror_type not in entries:
            continue
        lines.extend([
            f'{mirror_type} Runs',
            '-' * (len(mirror_type) + 5),
            '',
            f'.. toctree::',
            f'   :maxdepth: 1',
            '',
        ])
        for ts in entries[mirror_type]:
            lines.append(f'   {mirror_type}/{ts}')
        lines.append('')

    if not entries:
        lines.extend([
            'No run reports found yet. Run ``plot_ETM.py`` or ``plot_ITM.py``',
            'after an optimization to generate reports.',
            '',
        ])

    index_path = os.path.join(runs_dir, 'index.rst')
    with open(index_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def build_html():
    """Run ``make html`` in the docs/ directory to rebuild Sphinx pages.

    Failures are printed but do not raise — report .rst is still useful
    even if the HTML build fails (e.g. missing Sphinx install).
    """
    docs_dir = os.path.join(_PROJECT_ROOT, 'docs')
    if not os.path.isfile(os.path.join(docs_dir, 'Makefile')):
        return
    try:
        result = subprocess.run(
            ['make', 'html'],
            cwd=docs_dir,
            capture_output=True, text=True, timeout=120,
        )
        if result.returncode == 0:
            print(f"Sphinx HTML built: {docs_dir}/_build/html/")
        else:
            print(f"Sphinx build warning (exit {result.returncode}): "
                  f"{result.stderr[-200:] if result.stderr else 'no details'}")
    except FileNotFoundError:
        print("Sphinx build skipped: 'make' not found")
    except subprocess.TimeoutExpired:
        print("Sphinx build skipped: timed out after 120s")
