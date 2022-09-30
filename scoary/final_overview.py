import logging
from shutil import copy
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial.distance import cdist, squareform
import matplotlib as mpl

mpl.use('SVG')
# The SVG backend avoids this error message:
# ValueError: Image size of 700x165660 pixels is too large. It must be less than 2^16 in each direction.
# This allows for dendrograms with at least 20'000 traits

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.collections import PatchCollection, QuadMesh
from matplotlib.patches import Rectangle

from .utils import AnalyzeTraitNamespace, ROOT_DIR, RecursionLimit

logger = logging.getLogger('scoary.final_overview')


def plot_dendrogram(linkage_matrix: np.ndarray, labels: [str], ax: Axes) -> {}:
    with RecursionLimit(max(1000, len(linkage_matrix))):  # empirically tested for up to 20'000 traits
        dendrogram_params = hierarchy.dendrogram(
            linkage_matrix,
            orientation='left',
            labels=labels,
            no_labels=True,
            ax=ax
        )

    ax.set_xlim(left=1, right=0)
    ax.tick_params(
        axis='both', which='both',
        bottom=True, top=False, left=False, right=False,
        labelleft=False,
    )
    return dendrogram_params


def add_clickable_patches(patch_names, fig: Figure, ax: Axes, max_x: int):
    patches = [
        Rectangle(xy=(0, i * 10), width=max_x, height=10)
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


def plot_manhattan_like(qvals: pd.DataFrame, fig: Figure, ax: Axes, column_defs: {str: dict}) -> [QuadMesh]:
    # determine y intervals: [0, 10, 20, ...]
    y = np.arange(start=5, stop=len(qvals.index) * 10, step=10)

    ax.set_xlim(left=0)

    max_x = 1.
    for i, (col, def_) in enumerate(column_defs.items()):
        if col not in qvals.columns:
            continue

        # create pcolormesh with logarithmic scale
        x = -np.log10(qvals[col])

        max_x = max(max_x, x.max())

        ax.scatter(
            x, y,
            marker=def_['marker'],
            color=def_['color']
        )

    # add y ticks and labels
    ax.yaxis.tick_right()
    ytick_locations = np.arange(start=5, stop=len(qvals) * 10, step=10)
    ax.set_yticks(ytick_locations, qvals.index)
    ax.tick_params(
        axis='both', which='both',
        bottom=False, top=False, left=False, right=False,
        labelbottom=False
    )

    # add grid
    ax.set_axisbelow(True)
    ax.set_xticks(ticks=np.arange(0, max_x + 1, 1), minor=True)
    ax.grid(visible=True, which='both', axis='x', linestyle='dashed')

    # add shape on top of colormesh and ticks that can be made clickable
    add_clickable_patches(qvals.index, fig, ax, max_x)


def create_final_overview(summary_df: pd.DataFrame, ns: AnalyzeTraitNamespace, isolate_info_df: pd.DataFrame = None):
    logger.debug('Adding preliminary summary.tsv...')
    summary_df.to_csv(f'{ns.outdir}/summary.tsv', sep='\t')

    logger.info('Adding isolate info...')
    if isolate_info_df is not None:
        isolate_info_df.to_csv(f'{ns.outdir}/isolate_info.tsv', sep='\t')

    logger.info('Copying files...')
    # copy files
    for file in ['overview.html', 'trait.html']:
        copy(src=f'{ROOT_DIR}/templates/{file}', dst=f'{ns.outdir}/{file}')
    for file in ['config.json', 'trait.js', 'trait.css', 'overview.js', 'overview.css', 'favicon.svg']:
        copy(src=f'{ROOT_DIR}/templates/{file}', dst=f'{ns.outdir}/app/{file}')

    summary_df_index = list(summary_df.index)
    if len(summary_df) > 1:
        logger.info('Calculating dendrogram linkage matrix...')
        # prepare data, create linkage_matrix
        pre_jaccard = ns.traits_df[summary_df.index].astype('float').T
        pre_jaccard = ((pre_jaccard.fillna(0.5) * 2) - 1).astype('int')  # False -> -1, NAN -> 0, True -> 1

        # whether class=0 or class=1 is arbitrary. Calculate both possibilities, take minimum
        d1 = cdist(pre_jaccard, pre_jaccard, metric='jaccard')
        d2 = cdist(pre_jaccard, 0 - pre_jaccard, metric='jaccard')
        jaccard_distance = np.minimum(d1, d2) * 2  # multiply by 2 to make maximal distance 1 again

        jaccard_df = pd.DataFrame(jaccard_distance, columns=pre_jaccard.index, index=pre_jaccard.index)
        del d1, d2, jaccard_distance
        linkage_matrix = hierarchy.linkage(
            squareform(jaccard_df, checks=True),
            method='average',
            optimal_ordering=True
        )

        # calculate plot proportions
        content_height = max(3., len(jaccard_df) / 6)  # height dependent on number of compounds
        whitespace_abs = 0.6  # absolute amount of whitespace
        total_height = content_height + whitespace_abs  # total plot height
        whitespace_rel = whitespace_abs / 3 / total_height  # relative amount of whitespace

        # create matplotlib figure
        plt.close()
        fig = plt.figure(figsize=(8, total_height))  # , dpi=4)
        gs = fig.add_gridspec(
            nrows=1, ncols=2, width_ratios=(2, 1),
            left=0.05, right=0.6, bottom=whitespace_rel * 2, top=1 - whitespace_rel,
            wspace=0, hspace=0
        )

        # get axes objects with shared y-axis
        ax_dendrogram = fig.add_subplot(gs[0, 0])
        ax_colorbar = fig.add_subplot(gs[0, 1], sharey=ax_dendrogram)

        logger.info('Plotting dendrogram...')
        # plot dendrogram
        dendrogram_params = plot_dendrogram(linkage_matrix, labels=jaccard_df.columns.values, ax=ax_dendrogram)

        # reindex summary_df according to order in dendrogram
        summary_df_index = dendrogram_params['ivl']
        summary_df = summary_df.reindex(summary_df_index)

        column_defs = {
            'best_fisher_q': {'marker': '$f$', 'color': 'tab:green'},
            'best_empirical_p': {'marker': '$e$', 'color': 'tab:purple'},
            'best_fq*ep': {'marker': '*', 'color': 'tab:red'}
        }
        cols = [col for col in column_defs.keys() if col in summary_df.columns]
        plot_manhattan_like(summary_df[cols], fig=fig, ax=ax_colorbar, column_defs=column_defs)

        # save plot
        plt.savefig(f'{ns.outdir}/overview_plot.svg', format='svg')
        plt.close()

    if ns.trait_info_df is not None:
        logger.info('Adding trait info...')
        summary_df = summary_df \
            .merge(ns.trait_info_df, left_index=True, right_index=True, how='left', copy=False) \
            .reindex(summary_df_index)  # merging destroys index order

    # save summary_df, ensure order matches plot
    logger.info('Adding summary.tsv...')
    summary_df.index.name = 'Trait'
    summary_df.to_csv(f'{ns.outdir}/summary.tsv', sep='\t')
