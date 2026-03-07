"""Sphinx run report generation utilities.

Generates .rst files for individual optimization runs with inline plots,
summary tables, and quality assessment. Also maintains an auto-generated
index of all run reports.
"""

import os
import re
import shutil
import subprocess
import tempfile
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


def _is_ar(target):
    """Return True if target indicates anti-reflection (AR) coating."""
    return target > 0.5


def _format_trans_value(T_achieved, target):
    """Format a transmission metric value using AR/HR logic.

    AR targets (target > 0.5): show reflectivity in ppm.
    HR targets (target <= 0.5): show transmission in ppm.
    """
    if _is_ar(target):
        R = 1 - T_achieved
        return f'{R*1e6:.1f} ppm'
    return f'{T_achieved*1e6:.2f} ppm'


def _format_trans_target(target):
    """Format a transmission target using AR/HR logic."""
    if _is_ar(target):
        R_target = 1 - target
        return f'{R_target*1e6:.0f} ppm'
    return f'{target*1e6:.1f} ppm'


def _trans_label(cost_name, target, wavelength_nm):
    """Build a display label like 'T @ 1064 nm' or 'R @ 700 nm'."""
    prefix = 'R' if _is_ar(target) else 'T'
    return f'{prefix} @ {wavelength_nm:.0f} nm'


def _format_value(val, metric_name, target=None):
    """Format a metric value for display in the summary table."""
    if metric_name.startswith('Trans') and target is not None:
        return _format_trans_value(val, target)
    elif metric_name == 'Brownian':
        return f'{val:.4f}'
    elif metric_name == 'Thermooptic':
        return f'{val:.4e}'
    return f'{val:.4g}'


def _format_target(target, metric_name):
    """Format a target value for display."""
    if metric_name.startswith('Trans'):
        return _format_trans_target(target)
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
        for key in ('T1', 'T2', 'T3', 'scalarCost'):
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

        # Optimizer metadata
        try:
            metrics['algorithm'] = f['algorithm'][()].decode() \
                if isinstance(f['algorithm'][()], bytes) \
                else str(f['algorithm'][()])
        except (KeyError, TypeError):
            pass
        try:
            metrics['wall_time'] = float(np.array(f['wall_time']))
        except (KeyError, TypeError):
            pass

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
        'T1': 'Trans1',
        'T2': 'Trans2',
        'T3': 'Trans3',
    }

    # Compute actual wavelength labels from misc ratios
    misc = params.get('misc', {})
    wavelength = params.get('_wavelength', 1064e-9)
    wl_ratios = {'Trans1': 1.0,
                 'Trans2': misc.get('lambda2', 0.5),
                 'Trans3': misc.get('lambda3', None)}

    for hdf_key, cost_name in metric_to_cost.items():
        if hdf_key in metrics and cost_name in params.get('costs', {}):
            achieved = metrics[hdf_key]
            target = params['costs'][cost_name]['target']
            weight = params['costs'][cost_name]['weight']
            if weight == 0:
                continue
            status = assess_quality(achieved, target, cost_name)
            wl_ratio = wl_ratios.get(cost_name)
            if wl_ratio is not None:
                wl_nm = wavelength * wl_ratio * 1e9
                label = _trans_label(cost_name, target, wl_nm)
            else:
                label = cost_name
            lines.extend([
                f'   * - {label}',
                f'     - {_format_value(achieved, cost_name, target)}',
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
    if 'algorithm' in metrics:
        lines.extend([
            f'   * - Optimizer',
            f'     - {metrics["algorithm"]}',
            f'     -',
            f'     -',
        ])
    if 'wall_time' in metrics:
        lines.extend([
            f'   * - Wall time',
            f'     - {metrics["wall_time"]:.1f} s',
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

    if 'noise' in fig_paths:
        lines.extend([
            'Thermal Noise Budget',
            '--------------------', '',
            f'.. image:: {rel_prefix}/{fig_paths["noise"]}',
            '   :width: 100%', '',
        ])

    # MC section — build wavelength-aware labels dynamically
    # Row 0: T1 (PSL) always in ppm.  Row 1: T2.  Row 2 (optional): T3.
    # Use AR/HR logic: AR targets (>0.5) show R, HR targets show T.
    wl1_nm = wavelength * 1e9
    lambda2_ratio = misc.get('lambda2', 0.5)
    wl2_nm = wavelength * lambda2_ratio * 1e9
    lambda3_ratio = misc.get('lambda3', None)

    costs_cfg = params.get('costs', {})
    t2_target = costs_cfg.get('Trans2', {}).get('target', 0)
    t3_target = costs_cfg.get('Trans3', {}).get('target', 0)

    def _mc_trans_label(wl_nm, target):
        prefix = 'R' if _is_ar(target) else 'T'
        return rf':math:`{prefix}_{{{wl_nm:.0f}}}` [ppm]'

    # MC stores T [ppm] for HR and R [ppm] for AR — no transforms needed.
    mc_labels = [
        ('T1', rf':math:`T_{{{wl1_nm:.0f}}}` [ppm]'),
        ('T2', _mc_trans_label(wl2_nm, t2_target)),
    ]
    if lambda3_ratio is not None:
        wl3_nm = wavelength * lambda3_ratio * 1e9
        mc_labels.append(('T3', _mc_trans_label(wl3_nm, t3_target)))
    mc_labels.extend([
        ('S_TO', r':math:`S_\mathrm{TO}` '
         r'[:math:`\times 10^{-21}` m/:math:`\sqrt{\mathrm{Hz}}`]'),
        ('S_Br', r':math:`S_\mathrm{Br}` '
         r'[:math:`\times 10^{-21}` m/:math:`\sqrt{\mathrm{Hz}}`]'),
        ('E_surf', r':math:`E_\mathrm{surface}` [V/m]'),
    ])

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
    for entry in sorted(os.listdir(runs_dir)):
        subdir = os.path.join(runs_dir, entry)
        if not os.path.isdir(subdir):
            continue
        timestamps = []
        for fname in sorted(os.listdir(subdir), reverse=True):
            if fname.endswith('.rst') and fname != 'index.rst':
                timestamps.append(fname[:-4])
        if timestamps:
            entries[entry] = timestamps

    lines = [
        'Run Reports', '===========', '',
        'Auto-generated index of optimization run reports.', '',
    ]

    for mirror_type in sorted(entries):
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


def build_pdf(rst_path, output_path=None):
    """Build a single-page PDF from one run report RST.

    Creates a temporary Sphinx project containing just the target RST,
    copies referenced images so paths resolve, builds LaTeX, and
    compiles to PDF.

    Parameters
    ----------
    rst_path : str or Path
        Path to the run report ``.rst`` file (e.g.
        ``docs/runs/SFG/260306_021817.rst``).
    output_path : str or Path, optional
        Where to write the final PDF.  Defaults to the same directory
        as the RST with a ``.pdf`` extension.

    Returns
    -------
    str
        Path to the generated PDF.
    """
    rst_path = Path(rst_path).resolve()
    if not rst_path.exists():
        raise FileNotFoundError(f"RST not found: {rst_path}")

    if output_path is None:
        output_path = rst_path.with_suffix('.pdf')
    output_path = Path(output_path).resolve()

    # Read the RST and rewrite image paths to be relative to tmpdir
    rst_text = rst_path.read_text()

    tmpdir = Path(tempfile.mkdtemp(prefix='ob_pdf_'))
    try:
        # Find all image directives, convert SVG→PDF, rewrite paths
        img_pattern = re.compile(r'\.\. image:: (.+)')
        for m in img_pattern.finditer(rst_text):
            img_rel = m.group(1)
            img_src = (rst_path.parent / img_rel).resolve()
            if not img_src.exists():
                continue
            if img_src.suffix == '.svg':
                pdf_name = img_src.stem + '.pdf'
                img_dst = tmpdir / pdf_name
                subprocess.run(
                    ['rsvg-convert', '-f', 'pdf', '-o',
                     str(img_dst), str(img_src)],
                    check=True, capture_output=True,
                )
                rst_text = rst_text.replace(img_rel, pdf_name)
            else:
                img_dst = tmpdir / img_src.name
                shutil.copy2(img_src, img_dst)
                rst_text = rst_text.replace(img_rel, img_src.name)

        # Write the modified RST as index.rst
        (tmpdir / 'index.rst').write_text(rst_text)

        # Minimal conf.py for LaTeX/PDF build
        title = rst_path.stem
        conf = f"""\
project = 'OptimalBragg'
author = 'Caltech Experimental Gravity Group'
extensions = []
exclude_patterns = ['_build']
latex_elements = {{
    'papersize': 'letterpaper',
    'pointsize': '11pt',
}}
latex_documents = [
    ('index', 'report.tex', '{title}', author, 'howto'),
]
"""
        (tmpdir / 'conf.py').write_text(conf)

        # sphinx-build → LaTeX
        latex_dir = tmpdir / '_build' / 'latex'
        result = subprocess.run(
            ['sphinx-build', '-b', 'latex', str(tmpdir), str(latex_dir)],
            capture_output=True, text=True, timeout=120,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"Sphinx LaTeX build failed:\n{result.stderr[-500:]}")

        # Copy converted images into the latex build dir so pdflatex finds them
        for pdf_img in tmpdir.glob('*.pdf'):
            shutil.copy2(pdf_img, latex_dir / pdf_img.name)

        # latexmk → PDF
        result = subprocess.run(
            ['latexmk', '-pdf', '-interaction=nonstopmode', 'report.tex'],
            cwd=str(latex_dir),
            capture_output=True, text=True, timeout=120,
        )
        pdf_built = latex_dir / 'report.pdf'
        if not pdf_built.exists():
            raise RuntimeError(
                f"PDF compilation failed:\n{result.stderr[-500:]}")

        shutil.copy2(pdf_built, output_path)
        print(f"PDF saved: {output_path}")
        return str(output_path)

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
