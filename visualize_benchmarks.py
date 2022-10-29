import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

benchmarks_pd = pd.read_csv('benchmarks.csv', sep=';')

time_col = 'Elapsed Time'
iters_without_err_col = 'Iters without err counting'


def plot_benchmarks_by_target(target_device: str):
    versions = ['Without BLAS', 'BLAS naive',
                'BLAS without error copying', 'BLAS device pointer mode']
    benchmarks_by_device = benchmarks_pd[benchmarks_pd['Target device']
                                         == target_device]

    n_blocks_unique = np.unique(benchmarks_pd['Grid size'])

    fig, axs = plt.subplots(ncols=len(n_blocks_unique), figsize=(18, 9))

    axs[0].set_ylabel('Elapsed time')
    for ax_idx, n_blocks in enumerate(n_blocks_unique):
        ax = axs[ax_idx]
        benchmarks_by_blocks = benchmarks_by_device[benchmarks_pd['Grid size'] == n_blocks]

        bars = ax.bar(benchmarks_by_blocks['Algo ver'],
               benchmarks_by_blocks['Elapsed Time'])
        ax.set_xlabel(f'Grid size {n_blocks}')
        ax.set_xticks(list(range(len(benchmarks_by_blocks['Algo ver']))))
        ax.set_xticklabels(benchmarks_by_blocks['Algo ver'], rotation=10)
        ax.bar_label(bars)

    fig.suptitle(target_device)
    fig.savefig(f'benchmarks_{target_device}.png')


plot_benchmarks_by_target('CPU')
plot_benchmarks_by_target('GPU')
