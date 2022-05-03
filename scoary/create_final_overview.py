import json
from shutil import copy
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial.distance import cdist, squareform

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.collections import PatchCollection, QuadMesh
from matplotlib.patches import Rectangle
from matplotlib.colors import LogNorm, ListedColormap, Colormap, LinearSegmentedColormap

from .utils import MockNamespace, ROOT_DIR


def plot_dendrogram(linkage_matrix: np.ndarray, labels: [str], ax: Axes) -> {}:
    # no recursion problem here. Works with up to 20'000 isolates, becomes more of a memory problem.
    dendrogram_params = hierarchy.dendrogram(
        linkage_matrix,
        orientation='left',
        labels=labels,
        no_labels=True,
        ax=ax
    )

    ax.set_xlim([1, 0])
    ax.tick_params(
        axis='both', which='both',
        bottom=True, top=False, left=False, right=False,
        labelleft=False,
    )
    return dendrogram_params


def add_clickable_patches(patch_names, fig: Figure, ax: Axes):
    patches = [
        Rectangle(xy=(0, i * 10), width=15, height=10)
        for i in range(len(patch_names))
    ]

    # Create patch collection with specified colour/alpha
    pc = PatchCollection(
        patches,
        facecolors='white', alpha=1e-6,  # will result in 'opacity: 0'
        gid='clickable-patches',
        transform=ax.transData, figure=fig
    )

    # add urls -> this doesn't work
    # pc.set_urls([f'overview.html?trait={n}' for n in patch_names])

    fig.add_artist(pc)


def plot_qvals(qvals: pd.DataFrame, fig: Figure, ax: Axes, cmaps: {str: str | Colormap}) -> [QuadMesh]:
    # determine y intervals: [0, 10, 20, ...]
    y = np.arange(start=0, stop=len(qvals.index) * 10 + 1, step=10)

    pcms = []
    for i, (col, cmap) in enumerate(cmaps.items()):
        # determine x intervals: [0, 1] / [1, 2] / ...
        x = np.array([i, i + 1])

        # create pcolormesh with logarithmic scale
        pcm = ax.pcolormesh(
            x, y, qvals[[col]],
            cmap=cmap,
            norm=LogNorm(vmin=qvals[col].min(), vmax=1.)
        )
        pcms.append(pcm)

    # add ticks
    ytick_locations = np.arange(start=5, stop=len(qvals) * 10, step=10)
    ax.set_yticks(ytick_locations, qvals.index)
    ax.yaxis.tick_right()
    ax.tick_params(
        axis='both', which='both',
        bottom=False, top=False, left=False, right=False,
        labelbottom=False
    )

    # add shape on top of colormesh and ticks that can be made clickable
    add_clickable_patches(qvals.index, fig, ax)

    # return pcolormesh object, can be used to plot colorbar
    return pcms


def save_colorbars(pcms: [QuadMesh], cols: [str], out: str = None):
    fig = plt.figure(figsize=(len(pcms), 4))
    gs = fig.add_gridspec(
        nrows=1, ncols=len(pcms),
        left=0.05, right=0.99 - (0.5 / len(pcms)),
        bottom=0.01, top=0.93,
        wspace=1.5,
    )

    plt.rcParams['axes.titley'] = 1.05
    for i, pcm in enumerate(pcms):
        # add title on separate axis
        ax_title = fig.add_subplot(gs[i])
        ax_title.grid(False)
        plt.axis('off')
        ax_title.text(x=0.8, y=1.02, s=cols[i], fontsize=8, fontweight='bold', ha='center')

        # add colorbar
        ax_cbar = fig.add_subplot(gs[i])
        fig.colorbar(pcm, ax=None, cax=ax_cbar, extend='max')

    if out is None:
        plt.show()
    else:
        plt.savefig(out, format='svg')  # , bbox_inches='tight')
    plt.close()


def create_final_overview(overview_ds: pd.DataFrame, ns: MockNamespace, isolate_info_df: pd.DataFrame = None):
    # add isolate info
    if isolate_info_df is not None:
        isolate_info_df.to_csv(f'{ns.outdir}/isolate_info.tsv', sep='\t')

    # prepare data, create linkage_matrix
    pre_jaccard = ns.traits_df[overview_ds.index].astype('float').T
    pre_jaccard = ((pre_jaccard.fillna(0.5) * 2) - 1).astype('int')  # False -> -1, NAN -> 0, True -> 1

    # whether class=0 and class=1 are arbitrary. Calculate both possibilities, take minimum
    d1 = cdist(pre_jaccard, pre_jaccard, metric='jaccard')
    d2 = cdist(pre_jaccard, 0 - pre_jaccard, metric='jaccard')
    jaccard_distance = np.minimum(d1, d2)

    jaccard_df = pd.DataFrame(
        jaccard_distance, columns=pre_jaccard.index, index=pre_jaccard.index
    )
    linkage_matrix = hierarchy.linkage(squareform(jaccard_df), 'single')

    # calculate plot proportions
    content_height = max(3., len(jaccard_df) / 6)  # height dependent on number of compounds
    whitespace_abs = 0.6  # absolute amount of whitespace
    total_height = content_height + whitespace_abs  # total plot height
    whitespace_rel = whitespace_abs / 3 / total_height  # relative amount of whitespace

    # create matplotlib figure
    plt.close()
    fig = plt.figure(figsize=(7, total_height))
    gs = fig.add_gridspec(
        nrows=1, ncols=2, width_ratios=(7, 1),
        left=0.05, right=0.6, bottom=whitespace_rel * 2, top=1 - whitespace_rel,
        wspace=0, hspace=0,
        # wspace=0.05, hspace=0.05,
    )

    # get axes objects with shared y-axis
    ax_dendrogram = fig.add_subplot(gs[0, 0])
    ax_colorbar = fig.add_subplot(gs[0, 1], sharey=ax_dendrogram)

    # plot dendrogram
    dendrogram_params = plot_dendrogram(linkage_matrix, labels=jaccard_df.columns.values, ax=ax_dendrogram)

    # reindex overview_ds according to order in dendrogram
    dendrogram_index = dendrogram_params['ivl'][::-1]
    overview_ds = overview_ds.reindex(dendrogram_index)
    overview_ds.index.name = 'Trait'

    # plot qvals
    cmaps = {
        'min_qval': 'Spectral',
        'min_pval_empirical': LinearSegmentedColormap.from_list(
            name='pval_emp_cbar',
            colors=['#590d22', '#800f2f', '#a4133c', '#c9184a', '#ff4d6d',
                    '#ff758f', '#ff8fa3', '#ffb3c1', '#ffccd5', '#fff0f3'])
    }
    cols = [col for col in cmaps.keys() if col in overview_ds.columns]
    pcms = plot_qvals(overview_ds[cols], fig=fig, ax=ax_colorbar, cmaps=cmaps)

    # save plot
    plt.savefig(f'{ns.outdir}/overview_plot.svg', format='svg')
    plt.close()

    files = [f'{file_name}.{file_type}' for file_type in ('html', 'css', 'js') for file_name in ('overview', 'trait')]
    files.append('config.json')
    import os
    for file in files:
        os.symlink(src=f'{ROOT_DIR}/templates/{file}', dst=f'{ns.outdir}/{file}')
    copy(src=f'{ROOT_DIR}/templates/favicon.ico', dst=f'{ns.outdir}/favicon.ico')

    # create color bar, save
    save_colorbars(pcms, [c.removeprefix('min_') for c in cols], out=f'{ns.outdir}/overview_colorbar.svg')

    if ns.trait_info_df is not None:
        overview_ds = overview_ds \
            .merge(ns.trait_info_df, left_index=True, right_index=True, how='left', copy=False) \
            .reindex(dendrogram_index)  # merging destroys index order

    # save overview_ds (reversed index, so the order matches the plot)
    overview_ds.reindex(dendrogram_index).to_csv(f'{ns.outdir}/overview.tsv', sep='\t')
