from .progressbar import print_progress
from .utils import *
from .ScoaryTree import ScoaryTree
from .load_genes import load_genes
from .load_traits import load_traits
from .final_overview import create_final_overview
from .analyze_trait import analyze_trait, worker

logger = logging.getLogger('scoary')


def scoary(
        genes: str,
        traits: str,
        outdir: str,
        multiple_testing: str = 'bonferroni:0.999',
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
        random_state: int = None,
        limit_traits: (int, int) = None,
        version: bool = False  # Dummy variable, only used to create docstring (see main function)
) -> None:
    """
    Scoary2: Associate genes with traits!

    :param genes: Path to gene presence/absence table: columns=isolates, rows=genes
    :param traits: Path to trait presence/absence table: columns=traits, rows=isolates
    :param outdir: Directory to place output files
    :param multiple_testing: "method:cutoff" for filtering genes after Fisher's test, where cutoff is a number that
    specifies the FWER and method is one of [native, bonferroni, sidak, holm-sidak, holm, simes-hochberg, hommel,
    fdr_bh, fdr_by,  fdr_tsbh, fdr_tsbky]. Alternatively, the method can be 'native': then, the cutoff targets the
    uncorrected p-value from Fisher's test.
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
    :param random_state: Set a fixed seed for the random number generator
    :param limit_traits: Limit the analysis to traits n to m. Useful for debugging. Example: "(0, 10)"
    :param version: Print software version of Scoary2 and exit.
    """
    print(f'Welcome to Scoary2! ({get_version()})')

    # parse input, create outdir, setup logging
    trait_data_type = decode_unicode(trait_data_type)
    gene_data_type = decode_unicode(gene_data_type)
    if n_cpus_binarization is None:
        n_cpus_binarization = 1 + n_cpus // 10
    outdir = setup_outdir(outdir, input=locals())
    setup_logging(logger, f'{outdir}/logs/scoary-2.log')
    logger.debug(f'Scoary2 Version: {get_version()}')
    mt_f_method, mt_f_cutoff = parse_correction(multiple_testing)
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
    tree.write_newick(f'{outdir}/tree.nwk')

    all_labels = set(tree.labels())

    duplication_df = create_duplication_df(traits_df)

    logger.info('Finalizing setup...')
    if n_cpus == 1:
        ns, counter, lock = AnalyzeTraitNamespace(), MockCounter(), MockLock()
    else:
        from .init_multiprocessing import init, mp
        mgr, ns, counter, lock = init()

    ns = AnalyzeTraitNamespace.create_namespace(ns, {
        'start_time': datetime.now(),
        'counter': counter,
        'lock': lock,
        'outdir': outdir,
        'genes_orig_df': genes_orig_df,
        'genes_bool_df': genes_bool_df,
        'gene_info_df': gene_info,
        'numeric_df': numeric_df,
        'traits_df': traits_df,
        'trait_info_df': trait_info,
        'duplication_df': duplication_df,
        'tree': tree,
        'all_labels': all_labels,
        'mt_f_method': mt_f_method,
        'mt_f_cutoff': mt_f_cutoff,
        'max_genes': max_genes,
        'worst_cutoff': worst_cutoff,
        'n_permut': n_permut,
        'random_state': random_state,
        'pairwise': pairwise,
    })

    traits = traits_df.columns.to_list()

    logger.info('Start analyzing traits...')
    if n_cpus == 1:
        picking_start = datetime.now()
        trait_to_result = {trait: analyze_trait(trait, ns) for trait in traits}
    else:
        mp.freeze_support()
        queue = mgr.JoinableQueue()
        trait_to_result = mgr.dict()
        [queue.put(trait) for trait in traits]
        procs = [mp.Process(target=worker, args=(queue, ns, trait_to_result, i)) for i in range(n_cpus)]
        picking_start = datetime.now()
        [p.start() for p in procs]
        [p.join() for p in procs]
    picking_end = datetime.now()

    print_progress(
        len(ns.traits_df.columns), len(ns.traits_df.columns),
        message='COMPLETE!', start_time=ns.start_time, message_width=25,
        end='\n'
    )
    logging.info(f'Picking took {picking_end - picking_start}')
    try:
        summary_df = create_summary_df(trait_to_result)
        create_final_overview(summary_df, ns, isolate_info)

        logger.info('Complete success!')

    except NoTraitsLeftException as e:
        logger.info(str(e))

    logger.debug(f'Took {datetime.now() - start_time}')

    print(CITATION)


def create_summary_df(trait_to_result: {str: [dict | str | None]}) -> pd.DataFrame | None:
    """
    Turn trait_to_result into a pandas.DataFrame. Example:

             best_fisher_p  best_fisher_q  best_empirical_p  best_fq*ep
    Trait_1       0.574066   4.384058e-01          0.035964    0.035964
    Trait_2       0.432940   2.667931e-01          0.133866    0.133866
    Trait_3       0.194418   7.981206e-08          0.020979    0.691309

    :param trait_to_result: dictionary where keys are trait names and values are either dict|str|None
    :return: pandas.DataFrame
    """
    # res may contain: float, str or None
    # float: smallest empirical pvalue
    # str: duplicated trait -> name of trait
    # None: no gene was significant

    trait_to_result = {t: trait_to_result[r] if type(r) is str else r for t, r in
                       trait_to_result.items()}  # restore duplicates
    trait_to_result = {t: r for t, r in trait_to_result.items() if r is not None}  # remove Nones

    if len(trait_to_result) == 0:
        raise NoTraitsLeftException('No traits left after filtering')

    summary_df = pd.DataFrame(trait_to_result).T
    summary_df = summary_df.infer_objects()  # harmonize dtypes

    logger.debug(f'Created summary_df:\n{summary_df}')

    return summary_df


def create_duplication_df(traits_df: pd.DataFrame) -> pd.Series:
    """
    Returns a pd.Series that maps duplicated traits to the first occurrence
    """
    hash_df = pd.DataFrame(index=traits_df.columns)
    hash_df['hash'] = traits_df.apply(lambda x: hash(tuple(x)), axis=0)
    hash_df['is_duplicated'] = hash_df['hash'].duplicated(keep=False)
    hash_df['use_cache'] = hash_df['hash'].duplicated(keep='first')
    lookup_df = hash_df[hash_df['is_duplicated'] & ~hash_df['use_cache']].sort_values(by='hash')
    duplication_df = hash_df[hash_df['use_cache']]
    duplication_df = duplication_df['hash'].apply(
        func=lambda h: lookup_df.iloc[lookup_df.hash.searchsorted(h)].name
    )
    return duplication_df


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
