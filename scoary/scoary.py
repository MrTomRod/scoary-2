from scoary.progressbar import print_progress
from .utils import *
from .ScoaryTree import ScoaryTree
from .load_genes import load_genes
from .load_traits import load_traits
from .create_final_overview import create_final_overview
from .analyze_trait import analyze_trait, worker

logger = logging.getLogger('scoary-main')


def scoary(
        genes: str,
        traits: str,
        outdir: str,
        multiple_testing_fisher: str = 'bonferroni:0.999',
        multiple_testing_picking: str = 'bonferroni:0.999',
        gene_info: str = None,
        trait_info: str = None,
        isolate_info: str = None,
        newicktree: str = None,
        no_pairwise: bool = False,  # ?
        n_permut: int = 0,
        restrict_to: str = None,
        ignore: str = None,
        threads: int = 1,
        trait_data_type: str = 'binary:,',
        gene_data_type: str = 'gene-count',
        random_state: int = None,
        limit_traits: (int, int) = None,
):
    trait_data_type = decode_unicode(trait_data_type)
    gene_data_type = decode_unicode(gene_data_type)
    outdir = setup_outdir(outdir, input=locals())
    setup_logging(f'{outdir}/scoary-2.log')

    mt_f_method, mt_f_cutoff = parse_correction(multiple_testing_fisher)
    mt_p_method, mt_p_cutoff = parse_correction(multiple_testing_picking)

    assert n_permut == 0 or n_permut >= 100, f'{n_permut=} must be at least 100.'

    # load traits data  (numeric_df may be None)
    numeric_df, traits_df = load_traits(
        traits=traits,
        trait_data_type=trait_data_type,
        restrict_to=restrict_to,
        ignore=ignore,
        threads=threads,
        random_state=random_state,
        outdir=outdir,
        limit_traits=limit_traits
    )

    # dynamically set recursion limit, should work for ~ 13'000 isolates
    _recursion_limit = max(1000, 100 + len(traits_df.index) ** 2)
    logger.warning(f'Setting recursion limit to {_recursion_limit}')
    sys.setrecursionlimit(_recursion_limit)

    # load traits info
    trait_info_df = load_info_file(
        logger=logger, info_file=trait_info, merge_col='Trait',
        expected_overlap_set=set(traits_df.columns), reference_file=traits
    ) if trait_info else None

    # load genes data
    genes_orig_df, genes_bool_df = load_genes(
        genes,
        gene_data_type=gene_data_type,
        restrict_to=traits_df.index,
    )

    # load genes info
    gene_info_df = load_info_file(
        logger=logger, info_file=gene_info, merge_col='Gene',
        expected_overlap_set=set(genes_bool_df.index), reference_file=genes
    ) if gene_info else None

    # load isolate info
    isolate_info_df = load_info_file(
        logger=logger, info_file=isolate_info, merge_col='Isolate',
        expected_overlap_set=set(genes_bool_df.columns), reference_file='placeholder'
    ) if isolate_info else None

    # create convenient dictionary: {gene: {label: bool})
    all_label_to_gene = get_all_label_to_gene(genes_bool_df)

    # load phylogeny
    if newicktree is None:
        logger.info('Generating phylogenetic tree from gene presence-absence-matrix...')
        tree = ScoaryTree.from_presence_absence(genes_bool_df)
    else:
        with open(newicktree) as f:
            tree = ScoaryTree.from_newick(f.read())
    tree.write_newick(f'{outdir}/tree.nwk')

    all_labels = set(tree.labels())

    duplication_df = create_duplication_df(traits_df)

    if threads == 1:
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
        'gene_info_df': gene_info_df,
        'numeric_df': numeric_df,
        'traits_df': traits_df,
        'trait_info_df': trait_info_df,
        'duplication_df': duplication_df,
        'tree': tree,
        'all_labels': all_labels,
        'mt_f_method': mt_f_method,
        'mt_f_cutoff': mt_f_cutoff,
        'mt_p_method': mt_p_method,
        'mt_p_cutoff': mt_p_cutoff,
        'n_permut': n_permut,
        'all_label_to_gene': all_label_to_gene,
        'random_state': random_state,
        'no_pairwise': no_pairwise,
    })

    traits = traits_df.columns.to_list()

    if threads == 1:
        trait_to_result = {trait: analyze_trait(trait, ns) for trait in traits}
    else:
        mp.freeze_support()
        queue = mgr.JoinableQueue()
        trait_to_result = mgr.dict()
        [queue.put(trait) for trait in traits]
        procs = [mp.Process(target=worker, args=(queue, ns, trait_to_result, i)) for i in range(threads)]
        [p.start() for p in procs]
        [p.join() for p in procs]

    print_progress(
        len(ns.traits_df.columns), len(ns.traits_df.columns),
        message='COMPLETE!', start_time=ns.start_time, message_width=25,
        end='\n'
    )
    summary_df = create_summary_df(trait_to_result)
    if summary_df is None and len(summary_df) == 0:
        logging.warning(f'no content in summary_df')
    else:
        create_final_overview(summary_df, ns, isolate_info_df)

    print(CITATION)


def create_summary_df(trait_to_result: {str: [dict | str | None]}) -> pd.DataFrame | None:
    """
    Turn trait_to_result into a pandas.DataFrame
    :param trait_to_result:
    :return:
    """
    # res may contain: float, str or None
    # float: smallest empirical pvalue
    # str: duplicated trait -> name of trait
    # None: no gene was significant

    trait_to_result = {t: trait_to_result[r] if type(r) is str else r for t, r in
                       trait_to_result.items()}  # restore duplicates
    trait_to_result = {t: r for t, r in trait_to_result.items() if r is not None}  # remove Nones

    if len(trait_to_result) == 0:
        logging.warning(f'No traits left after filtering!')
        return

    summary_df = pd.DataFrame(trait_to_result).T
    """      min_pval      min_qval  min_pval_empirical  min_qval_empirical
    Trait_1  0.574066  4.384058e-01            0.035964            0.035964
    Trait_2  0.432940  2.667931e-01            0.133866            0.133866
    Trait_3  0.194418  7.981206e-08            0.020979            0.691309
    ...
    """

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


CITATION = '''
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


If you use Scoary 2, please cite:
Roder T, Brynildsrud O, ... Title title title title title title
title title title title title title.
Journal. 2022;00:000.
'''.strip('\n')


def main():
    import fire

    fire.Fire(scoary)


if __name__ == '__main__':
    main()
