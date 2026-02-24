"""Background MC pipeline: mc → corner → rebuild report.

Runs as a subprocess launched by ``optimalbragg run``::

    python -m OptimalBragg.mc_pipeline <layers.hdf5> <mc_output.hdf5> <n_samples> <params.yml>
"""

import sys


def main():
    if len(sys.argv) < 5:
        print("Usage: python -m OptimalBragg.mc_pipeline "
              "<layers.hdf5> <mc_output.hdf5> <n_samples> <params.yml>")
        sys.exit(1)

    layers_hdf5 = sys.argv[1]
    mc_output = sys.argv[2]
    n_samples = int(sys.argv[3])
    params_yml = sys.argv[4]

    # 1. Run MC
    from OptimalBragg.mc import run_mc, save_mc
    print(f"MC pipeline: running {n_samples} samples...")
    result = run_mc(layers_hdf5, n_samples=n_samples)
    save_mc(result, mc_output)
    print(f"MC pipeline: saved {mc_output}")

    # 2. Rebuild plots + report with MC data (includes corner plot)
    from OptimalBragg.cli import _generate_plots_and_report
    _generate_plots_and_report(params_yml, hdf5_path=layers_hdf5,
                               mc_hdf5_path=mc_output)
    print("MC pipeline: report rebuilt with MC stats + corner plot.")


if __name__ == '__main__':
    main()
