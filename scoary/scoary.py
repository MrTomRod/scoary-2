from .progressbar import print_progress
from .utils import *
from .ScoaryTree import ScoaryTree
from .load_genes import load_genes
from .load_traits import load_traits
from .final_overview import create_final_overview
from .analyze_trait import analyze_trait_step_1_fisher, analyze_trait_step_2_pairpicking, worker, multiple_testing_correction

logger = logging.getLogger('scoary')


def scoary(
        genes: str,
        traits: str,
        outdir: str,
        multiple_testing: str = 'bonferroni:0.999',
        trait_wise_correction: bool = False,
        worst_cutoff: float = None,
        max_genes: int = None,
        gene_info: str = None,
        trait_info: str = None,
        isolate_info: str = None,
        newicktree: str = None,
        pairwise: bool = True,
        n_permut: int = 500,
        restrict_to: str = None,
        ignore: str = None,
        n_cpus: int = 1,
        n_cpus_binarization: int = None,
        trait_data_type: str = 'binary:,',
        gene_data_type: str = 'gene-count:,',
        force_binary_clustering: bool = False,
        symmetric: bool = True,
        distance_metric: str = 'jaccard',
        linkage_method: str = 'ward',
        optimal_ordering: bool = True,
        corr_method: str = 'pearson',
        random_state: int = None,
        limit_traits: (int, int) = None,
        version: bool = False  # Dummy variable, only used to create docstring (see main function)
) -> None:
    """
    Scoary2: Associate genes with traits!

    :param genes: Path to gene presence/absence table: columns=isolates, rows=genes
    :param traits: Path to trait presence/absence table: columns=traits, rows=isolates
    :param outdir: Directory to place output files
    :param multiple_testing: Apply multiple testing to the p-values of Fisher's test to account for the many
    genes/traits tested. Format: "method:cutoff".
    Cutoff is a number that specifies the FWER and method is one of [native, bonferroni, sidak, holm-sidak, holm,
    simes-hochberg, hommel, fdr_bh, fdr_by,  fdr_tsbh, fdr_tsbky].
    If method is 'native': then, the cutoff targets the uncorrected p-value from Fisher's test.
    :param trait_wise_correction: Apply multiple testing correction to each trait separately. Not recommended as
    this can lead to many false positives!
    :param worst_cutoff: Drop traits if no gene with "worst" p-value lower than threshold. Recommended if
    dataset contains multiple species
    :param max_genes: Keep only n highest-scoring genes in Fisher's test. Recommended if dataset is big and contains
     multiple species; avoids waisting computational resources on traits that simply correlate with phylogeny
    :param gene_info: Path to file that describes genes: columns=arbitrary properties, rows=genes
    :param trait_info: Path to file that describes traits: columns=arbitrary properties, rows=traits
    :param isolate_info: Path to file that describes isolates: columns=arbitrary properties, rows=isolates
    :param newicktree: Path to a custom tree in Newick format
    :param pairwise: If False, only perform Fisher's test. If True, also perform pairwise comparisons
     algorithm.
    :param n_permut: Post-hoc label-switching test: perform N permutations of the phenotype by random label switching.
     Low p-values suggest that the effect is not merely lineage-specific.
    :param restrict_to: Comma-separated list of isolates to which to restrict this analysis
    :param ignore: Comma-separated list of isolates to be ignored for this analysis
    :param n_cpus: Number of CPUs that should be used. There is overhead in multiprocessing, so if the dataset is
    small, use n_cpus=1
    :param n_cpus_binarization: Number of CPUs that should be used for binarization. Default: one tenth of n_cpus
    :param trait_data_type: "<method>:<?cutoff>:<?covariance_type>:<?alternative>:<?delimiter>" How to read the traits
     table. Example: "gene-list:\\t" for OrthoFinder N0.tsv table
    :param gene_data_type: "<data_type>:<?delimiter>" How to read the genes table. Example: "gene-list:\\t" for
     OrthoFinder N0.tsv table
    :param force_binary_clustering: Force clustering of binary data even if numeric data is available
    :param symmetric: if True, correlated and anti-correlated traits will cluster together
    :param distance_metric: distance metric (binary data only); See metric in https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.cdist.html
    :param linkage_method: linkage method for clustering [single, complete, average, weighted, ward, centroid, median]
    :param optimal_ordering: whether to use optimal ordering; See scipy.cluster.hierarchy.linkage.
    :param corr_method: correlation method (numeric data only) [pearson, kendall, spearman]
    :param random_state: Set a fixed seed for the random number generator
    :param limit_traits: Limit the analysis to traits n to m. Useful for debugging. Example: "(0, 10)"
    :param version: Print software version of Scoary2 and exit.
    """
    SCOARY_PRINT_CITATION = os.environ.get('SCOARY_PRINT_CITATION', 'TRUE') == 'TRUE'
    if SCOARY_PRINT_CITATION:
        print(f'Welcome to Scoary2! ({get_version()})')

    # parse input, create outdir, setup logging
    trait_data_type = decode_unicode(trait_data_type)
    gene_data_type = decode_unicode(gene_data_type)
    if n_cpus_binarization is None:
        n_cpus_binarization = 1 + n_cpus // 10
    outdir = setup_outdir(outdir, input=locals())

    setup_logging(logger, f'{outdir}/logs/scoary-2.log')

    logger.debug(f'Scoary2 Version: {get_version()}')
    mt_f_method, mt_f_cutoff = parse_correction(multiple_testing, 'multiple_testing')
    assert n_permut == 0 or n_permut >= 100, f'{n_permut=} must be at least 100.'

    # start
    start_time = datetime.now()

    # load traits data  (numeric_df may be None)
    logger.info('Loading traits...')
    numeric_df, traits_df = load_traits(
        traits=traits,
        trait_data_type=trait_data_type,
        restrict_to=restrict_to,
        ignore=ignore,
        n_cpus=n_cpus_binarization,
        random_state=random_state,
        outdir=outdir,
        limit_traits=limit_traits
    )

    # dynamically set recursion limit, should work for ~ 13'000 isolates
    _recursion_limit = max(1000, 100 + len(traits_df.index) ** 2)
    logger.debug(f'Setting recursion limit to {_recursion_limit}')
    sys.setrecursionlimit(_recursion_limit)

    if trait_info:
        logger.info('Loading trait info...')
        trait_info = load_info_file(
            logger=logger, info_file=trait_info, merge_col='Trait',
            expected_overlap_set=set(traits_df.columns), reference_file=traits
        )

    logger.info('Loading genes...')
    genes_orig_df, genes_bool_df = load_genes(
        genes,
        gene_data_type=gene_data_type,
        restrict_to=traits_df.index,
    )

    if gene_info:
        logger.info('Loading gene info...')
        gene_info = load_info_file(
            logger=logger, info_file=gene_info, merge_col='Gene',
            expected_overlap_set=set(genes_bool_df.index), reference_file=genes
        )

    if isolate_info:
        logger.info('Loading isolate info...')
        isolate_info = load_info_file(
            logger=logger, info_file=isolate_info, merge_col='Isolate',
            expected_overlap_set=set(genes_bool_df.columns), reference_file='placeholder'
        )

    # load phylogeny
    if newicktree is None:
        logger.info('Generating phylogenetic tree from gene presence-absence-matrix...')
        tree = ScoaryTree.from_presence_absence(genes_bool_df)
    else:
        logger.info('Loading phylogenetic tree from newick file...')
        with open(newicktree) as f:
            tree = ScoaryTree.from_newick(f.read())
        tree = tree.prune(genes_bool_df.columns)
    tree.write_newick(f'{outdir}/tree.nwk')

    all_labels = set(tree.labels())

    traits = traits_df.columns.to_list()
    duplicates = find_duplicates(traits_df)

    logger.info('Finalizing setup...')
    if n_cpus == 1:
        ns, counter, lock = AnalyzeTraitNamespace(), MockCounter(), MockLock()
    else:
        from .init_multiprocessing import init, mp
        mgr, ns, counter, lock = init()

    ns = AnalyzeTraitNamespace.create_namespace(ns, {
        'start_time': datetime.now(),
        'counter': counter,
        'queue_size': len(traits),
        'lock': lock,
        'outdir': outdir,
        'genes_orig_df': genes_orig_df,
        'genes_bool_df': genes_bool_df,
        'gene_info_df': gene_info,
        'numeric_df': numeric_df,
        'traits_df': traits_df,
        'trait_info_df': trait_info,
        'duplicates': duplicates,
        'tree': tree,
        'all_labels': all_labels,
        'mt_f_method': mt_f_method,
        'mt_f_cutoff': mt_f_cutoff,
        'trait_wise_correction': trait_wise_correction,
        'max_genes': max_genes,
        'worst_cutoff': worst_cutoff,
        'n_permut': n_permut,
        'random_state': random_state,
        'pairwise': pairwise,
        'multiple_testing_df': None,
    })

    logger.info('Starting step 1: Fisher\'s test...')
    if n_cpus == 1:
        step_1_start = datetime.now()
        trait_to_result = {trait: analyze_trait_step_1_fisher(trait, ns) for trait in traits}
    else:
        mp.freeze_support()
        queue = mgr.JoinableQueue()
        trait_to_result = mgr.dict()
        [queue.put(trait) for trait in traits]
        procs = [mp.Process(target=worker, args=(queue, ns, 1, trait_to_result, i)) for i in range(n_cpus)]
        step_1_start = datetime.now()
        [p.start() for p in procs]
        [p.join() for p in procs]

    step_1_end = datetime.now()
    print_progress(
        len(traits), len(traits),
        message='Step 1 complete!', start_time=step_1_start, message_width=25,
        end='\n'
    )
    logger.info(f'Step 1 took {step_1_end - step_1_start}')

    duplicated_traits = {trait: res for trait, res in trait_to_result.items() if type(res) is str}
    logger.info(f'Number of duplicated traits: {len(duplicated_traits)}')
    logger.info(f'Number of non-duplicated traits: {len(trait_to_result) - len(duplicated_traits)}')

    # multiple testing correction
    if trait_wise_correction:
        traits_left = {trait for trait, res in trait_to_result.items() if res is True}
        ns.multiple_testing_df = 'Not used'
    else:
        trait_to_result = {trait: res for trait, res in trait_to_result.items() if type(res) is not str}
        multiple_testing_df = multiple_testing_correction(
            pd.concat(trait_to_result), 'fisher_p', 'fisher_q',
            ns.mt_f_method, ns.mt_f_cutoff, False
        )
        multiple_testing_df.drop('fisher_p', axis=1, inplace=True)
        traits_left = multiple_testing_df.index.get_level_values(0).unique().to_list()
        ns.multiple_testing_df = multiple_testing_df
    del trait_to_result

    # Step 2: Pairpicking
    ns.queue_size = len(traits_left)
    ns.counter.value = 0
    logger.info(f'Number of traits left after multiple testing correction: {len(traits_left)}')

    logger.info('Starting step 2: Pair picking...')
    if n_cpus == 1:
        step_2_start = datetime.now()
        trait_to_result = {trait: analyze_trait_step_2_pairpicking(trait, ns) for trait in traits_left}
    else:
        mp.freeze_support()
        queue = mgr.JoinableQueue()
        trait_to_result = mgr.dict()
        [queue.put(trait) for trait in traits_left]
        procs = [mp.Process(target=worker, args=(queue, ns, 2, trait_to_result, i)) for i in range(n_cpus)]
        step_2_start = datetime.now()
        [p.start() for p in procs]
        [p.join() for p in procs]

    step_2_end = datetime.now()
    print_progress(
        len(traits_left), len(traits_left),
        message='Step 2 complete!', start_time=step_2_start, message_width=25,
        end='\n'
    )
    logger.info(f'Step 2 took {step_2_end - step_2_start}')

    try:
        summary_df = create_summary_df(trait_to_result, duplicated_traits)
    except NoTraitsLeftException as e:
        logger.info(str(e))
        logger.debug(f'Took {datetime.now() - start_time}')
        return
    del trait_to_result

    summary_df = summary_df.sort_values(
        by='best_fq*ep' if 'best_fq*ep' in summary_df.columns else 'best_fisher_q',
        ascending=False
    )

    create_final_overview(summary_df, ns.traits_df, ns.numeric_df, ns.outdir, ns.trait_info_df, isolate_info,
                          force_binary_clustering, symmetric, distance_metric, linkage_method, optimal_ordering, corr_method)

    logger.info('Cleaning up...')
    clean_up(outdir, summary_df.index.to_list())

    logger.info('Complete success!')

    logger.info(f'Took {datetime.now() - start_time}')

    if SCOARY_PRINT_CITATION:
        print(CITATION)


def create_summary_df(trait_to_result: {str: [dict | None]}, duplicated_traits: {str: str}) -> pd.DataFrame | None:
    """
    Turn trait_to_result into a pandas.DataFrame. Example:

             best_fisher_p  best_fisher_q  best_empirical_p  best_fq*ep
    Trait_1       0.574066   4.384058e-01          0.035964    0.035964
    Trait_2       0.432940   2.667931e-01          0.133866    0.133866
    Trait_3       0.194418   7.981206e-08          0.020979    0.691309

    :param trait_to_result: dictionary where keys are trait names and values are either dict|str|None
    :return: pandas.DataFrame
    """
    # res may contain: dict or None:
    #  - dict:  data to be added to summary_df as a row
    #  - None:  no gene was significant

    # remove Nones
    trait_to_result = {t: r for t, r in trait_to_result.items() if r is not None}

    # remove traits with no significant genes
    trait_to_result.update({t: trait_to_result[r] for t, r in duplicated_traits.items() if r in trait_to_result})

    if len(trait_to_result) == 0:
        raise NoTraitsLeftException('No traits left after filtering')

    summary_df = pd.DataFrame(trait_to_result).T
    summary_df = summary_df.infer_objects()  # harmonize dtypes

    logger.debug(f'Created summary_df:\n{summary_df}')

    return summary_df


def find_duplicates(traits_df: pd.DataFrame) -> pd.Series:
    """
    Returns a pd.Series that maps duplicated traits to the first occurrence
    """
    hash_df = pd.DataFrame(index=traits_df.columns)
    hash_df['hash'] = traits_df.apply(lambda x: hash(tuple(x)), axis=0)
    hash_df['is_duplicated'] = hash_df['hash'].duplicated(keep=False)
    hash_df['use_cache'] = hash_df['hash'].duplicated(keep='first')
    lookup_df = hash_df[hash_df['is_duplicated'] & ~hash_df['use_cache']].sort_values(by='hash')
    duplicates = hash_df[hash_df['use_cache']]
    duplicates = duplicates['hash'].apply(
        func=lambda h: lookup_df.iloc[lookup_df.hash.searchsorted(h)].name
    )
    return duplicates


def clean_up(outdir: str, traits_left: list[str]) -> None:
    import shutil
    for trait in os.listdir(f'{outdir}/traits'):
        if trait not in traits_left:
            shutil.rmtree(f'{outdir}/traits/{trait}')


CITATION = f'''
  ██████  ▄████▄   ▒█████   ▄▄▄       ██▀███ ▓██   ██▓   ░▒█████▒░ 
▒██    ▒ ▒██▀ ▀█  ▒██▒  ██ ▒████▄    ▓██   ██ ▒██  ██▒   ▒█▒   ██▒░
░ ▓██▄   ▒▓█    ▄ ▒██░  ██ ▒██  ▀█▄  ▓██ ░▄█   ▒██ ██░       ░█▀   
  ▒   ██ ▒▓▓▄ ▄██▒▒██   ██ ░██▄▄▄▄██ ▒██▀▀█▄   ░ ▐██▓░      ▄█     
▒██████▒   ▓███▀ ░░ ████▓▒░ ▓█   ▓██▒░██▓ ▒██▒   ██▒▓░   ░███████▒ 
▒ ▒▓▒ ▒ ░  ░▒ ▒  ░░ ▒░▒░▒░  ▒▒   ▓▒█░░ ▒▓ ░▒▓░  ██▒▒▒    ░▒▒  ░▒░  
░ ░▒  ░ ░  ░  ▒     ░ ▒ ▒░   ▒   ▒▒ ░  ░▒ ░ ▒░▓██ ░▒░     ░░   ▒░  
░  ░  ░  ░        ░ ░ ░ ▒    ░   ▒     ░░   ░ ▒ ▒ ░░           ░   
      ░  ░ ░          ░ ░        ░  ░   ░     ░ ░                  
         ░                                    ░ ░                  
                        Microbial Pan-GWAS


If you use Scoary2 ({get_version()}), please cite:
Roder, T. et al. Scoary2: Rapid association of phenotypic multi-omics 
data with microbial pan-genomes.
BioRxiv (2023) doi:10.1101/2023.04.19.537353.
'''.strip('\n')


def main():
    import sys, fire

    if '--version' in sys.argv:
        print(f'{get_version()}')
        exit(0)

    fire.Fire(scoary)


if __name__ == '__main__':
    main()
