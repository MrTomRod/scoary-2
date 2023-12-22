import logging
import os.path

import pandas as pd
import matplotlib as mpl

mpl.use('SVG')
# The SVG backend avoids this error message:
# ValueError: Image size of 700x165660 pixels is too large. It must be less than 2^16 in each direction.
# This allows for dendrograms with at least 20'000 traits

import mgwas_data_exploration_app.main as exploration_app

logger = logging.getLogger('scoary.final_overview')

SCORES_CONFIG = {
    "best_fisher_q": {
        "legend": "Fisher's <i>q</i>-value",
        "marker-matplotlib": "$f$",
        "marker-html": "<i>f</i>",
        "color": "forestgreen"
    },
    "best_empirical_p": {
        "legend": "Empirical <i>p</i>-value",
        "marker-matplotlib": "$e$",
        "marker-html": "<i>e</i>",
        "color": "mediumpurple"
    },
    "best_fq*ep": {
        "legend": "<i>fq*ep</i> score",
        "marker-matplotlib": "*",
        "marker-html": "*",
        "color": "crimson"
    }
}


def create_final_overview(
        summary_df: pd.DataFrame,
        traits_df: pd.DataFrame,
        numeric_df: pd.DataFrame,
        outdir: str,
        trait_info_df: pd.DataFrame = None,
        isolate_info_df: pd.DataFrame = None,
        force_binary_clustering: bool = False,
        symmetric: bool = True,
        distance_metric: str = 'jaccard',
        linkage_method: str = 'ward',
        optimal_ordering: bool = True,
        corr_method: str = 'pearson'
):
    # copy files from exploration app
    logger.info('Copying exploration app...')
    exploration_app.copy_app(outdir, config={'scores': SCORES_CONFIG})

    if isolate_info_df is not None:
        logger.info('Adding isolate_info.tsv...')
        isolate_info_df.to_csv(f'{outdir}/isolate_info.tsv', sep='\t')

    logger.debug('Adding preliminary summary.tsv...')
    summary_df.index.name = 'Trait'
    summary_df.to_csv(f'{outdir}/summary_orig.tsv', sep='\t')

    # append trait info
    if trait_info_df is not None:
        logger.debug('Adding trait_info_df to summary.tsv...')
        summary_df_index = list(summary_df.index)
        summary_df = summary_df \
            .merge(trait_info_df, left_index=True, right_index=True, how='left', copy=False) \
            .reindex(summary_df_index)  # merging destroys index order
        summary_df.index.name = 'Trait'
        summary_df.to_csv(f'{outdir}/summary.tsv', sep='\t')

    if len(summary_df) > 1:
        logger.info('Calculating dendrogram linkage matrix...')
        if numeric_df is None or force_binary_clustering:
            logger.info(f'Calculating dendrogram based on binary data using jaccard distances...')
            linkage_matrix, labels = exploration_app.calculate_linkage_matrix_from_binary(
                summary_df=summary_df,
                traits_df=traits_df,
                symmetric=symmetric,
                distance_metric=distance_metric,
                linkage_method=linkage_method,
                optimal_ordering=optimal_ordering
            )
        else:
            logger.info(f'Calculating dendrogram based on correlation of numeric features...')
            linkage_matrix, labels = exploration_app.calculate_linkage_matrix_from_numeric(
                summary_df=summary_df,
                traits_df=numeric_df,
                symmetric=symmetric,
                scale=True,
                corr_method=corr_method,
                linkage_method=linkage_method,
                optimal_ordering=optimal_ordering,
            )

        logger.info('Calculating dendrogram plot...')
        summary_df = exploration_app.final_plot(
            linkage_matrix=linkage_matrix,
            labels=labels,
            summary_df=summary_df,
            scores_config=SCORES_CONFIG,
            workdir=outdir,
            dendrogram_x_scale='linear',
            scores_x_scale='manhattan'
        )

        # save summary_df, ensure order matches plot
        logger.info('Saving sorted summary.tsv...')
        summary_df.index.name = 'Trait'
        summary_df.to_csv(f'{outdir}/summary.tsv', sep='\t')

    if not os.path.isfile(f'{outdir}/summary.tsv'):
        logger.debug('Moving summary_orig.tsv to summary.tsv...')
        os.rename(f'{outdir}/summary_orig.tsv', f'{outdir}/summary.tsv')
