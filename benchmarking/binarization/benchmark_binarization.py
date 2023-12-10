import os
import numpy as np
import pandas as pd
from scipy.stats import norm
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns


def create_common_ancestor_genome(core_genes: int = 3000, pan_genes: int = 6000, causal_genes: [float] = [0.01]) -> np.array:
    return np.concatenate((
        np.full(core_genes, True, dtype=bool),  # core genes are always present
        np.full(pan_genes + len(causal_genes), False, dtype=bool)  # pan genes and causal genes are initially absent
    ))


def create_mut_chance_series(core_genes: int = 3000, pan_genes: int = 6000, causal_genes: [float] = [0.01]) -> np.array:
    return np.concatenate((
        np.zeros(core_genes),  # core genes have no mutation chance
        np.random.rand(pan_genes) / 100,  # each pan gene has a random mutation chance between 0 and 0.01
        np.array(causal_genes)  # causal genes have specified mutation chance
    )).reshape(core_genes + pan_genes + len(causal_genes), 1)  # create 2D array


def mutate_genomes(genomes: pd.DataFrame, mut_change_series: np.array):
    """Mutate genomes"""
    random_values = np.random.rand(*genomes.shape)
    mutated = random_values <= mut_change_series
    # genome_collection and mutated are arrays of the same size. Wherever mutated is True, flip the bit in genome_collection
    return pd.DataFrame(
        np.logical_xor(genomes.values, mutated),
        index=genomes.index,
        columns=genomes.columns
    )


def branch_genomes(genomes: pd.DataFrame, mut_change_series: np.array, branch_probability: float):
    """Branch genomes"""
    random_values = np.random.rand(len(genomes.columns))
    branch = random_values < branch_probability
    if branch.sum() == 0:  # no new genomes
        return genomes
    new_genomes = genomes.loc[:, branch].copy()  # .copy necessary?
    column_names = [f'genome_{i}' for i in range(len(genomes.columns), len(genomes.columns) + len(new_genomes.columns))]
    new_genomes.columns = column_names
    mutate_genomes(new_genomes, mut_change_series=mut_change_series)
    return pd.concat([genomes, new_genomes], axis=1)


def create_genomes(n_genomes: int, core_genes: int = 3000, pan_genes: int = 6000, causal_genes: [float] = [0.01]):
    mut_change_series = create_mut_chance_series(core_genes, pan_genes, causal_genes)

    genomes = pd.DataFrame(
        data={'genome_0': create_common_ancestor_genome(core_genes, pan_genes, causal_genes)},
        index=[f'core_{x:05}' for x in range(core_genes)] +
              [f'pan_{x:05}' for x in range(pan_genes)] +
              [f'causal_{x:05}' for x in range(len(causal_genes))],
        dtype="bool"
    )

    while len(genomes.columns) < n_genomes:
        # mutate all genomes
        genomes = mutate_genomes(genomes, mut_change_series)
        # create new genomes/branches
        genomes = branch_genomes(genomes, mut_change_series=mut_change_series, branch_probability=0.01)

    # The last iteration of branch_genomes may have created more genomes than necessary.
    # Return only first n_genomes genomes
    genomes = genomes.iloc[:, :n_genomes]
    return genomes


def calculate_phenotype(genomes: pd.DataFrame, effect_size: float):
    """
    Calculate phenotype asuming a normal distribution.
    """
    return pd.Series(
        np.random.normal(
            loc=genomes.loc['causal_00000'].astype(int) * effect_size,
            scale=1
        ),
        index=genomes.columns,
        name='phenotype'
    )


def write_files(genomes: pd.DataFrame, genomes_file: str, phenotype: pd.Series, phenotype_file: str):
    genomes.astype('int').to_csv(genomes_file, sep='\t')
    phenotype.to_csv(phenotype_file, sep='\t')


def test(n_genomes: int, effect_size: float, core_genes: int = 3000, pan_genes: int = 6000, causal_genes: [float] = [0.01]):
    genomes = create_genomes(n_genomes, core_genes, pan_genes, causal_genes)
    print(genomes)
    phenotype = calculate_phenotype(genomes, effect_size)
    print(phenotype)
    write_files(genomes, 'simulations/genomes.tsv', phenotype, 'simulations/phenotype.tsv')


def _simulate(n_replicates, n, e, r):
    if os.path.isdir(f'simulations/{n=}_{e=}_{r=}'):
        return

    print(f'Running replicate {r} of {n_replicates} for {n} genomes and {e} effect size')
    np.random.seed(n + int(e * 2) + r)  # dirty hack: each replicate gets a predictable seed

    genomes = create_genomes(n)
    phenotype = calculate_phenotype(genomes, e)

    os.makedirs(f'simulations/{n=}_{e=}_{r=}', exist_ok=True)
    write_files(
        genomes, f'simulations/{n=}_{e=}_{r=}/genomes.tsv',
        phenotype, f'simulations/{n=}_{e=}_{r=}/phenotype.tsv'
    )


def generate_simulations():
    n_replicates = 20
    n_genomes = [25, 50, 75, 100, 150, 200]
    effect_size = [0.5, 0.75, 1., 1.5, 2., 2.5, 3.]

    os.makedirs('simulations', exist_ok=True)

    from multiprocessing import Pool
    with Pool() as pool:
        pool.starmap(
            _simulate,
            [(n_replicates, n, e, r) for n in n_genomes for e in effect_size for r in range(n_replicates)]
        )


def _scoary(msg: str, simulation: str):
    os.environ['SCOARY_RESET_LOGGERS'] = 'TRUE'
    os.environ['SCOARY_LOGLEVEL_STDOUT'] = 'WARNING'
    os.environ['SCOARY_PRINT_CITATION'] = 'FALSE'
    os.environ['SCOARY_PRINT_PROGRESS'] = 'FALSE'

    from scoary import scoary

    print(f'{msg}: Analyzing {simulation}')

    genes = f'simulations/{simulation}/genomes.tsv'
    traits = f'simulations/{simulation}/phenotype.tsv'
    outdir = f'simulations/{simulation}/scoary'

    if os.path.isdir(outdir):
        return
        # import shutil
        # shutil.rmtree(outdir)

    for file in [genes, traits]:
        assert os.path.exists(file), f'{file} does not exist'

    scoary(
        genes,
        traits,
        outdir,
        trait_data_type='gaussian:kmeans:\t',
        gene_data_type='gene-count:\t',
        multiple_testing='native:0.05',
        n_permut=1000,
        n_cpus=1,
        random_state=42,
    )

    assert os.path.isdir(outdir), f'{outdir} does not exist'

    if not os.listdir(f'{outdir}/traits'):
        print(f'{simulation=}: No traits found')


def analyze_scoary_results():
    datapoints = []
    for simulation in os.listdir('simulations'):
        datapoint = {key: float(value) if '.' in value else int(value)
                     for key, value in (pair.split('=') for pair in simulation.split('_'))}
        try:
            df = pd.read_csv(f'simulations/{simulation}/scoary/traits/phenotype/result.tsv', sep='\t', index_col=0)
            assert 'causal_00000' in df.index, f'{simulation=}: causal_00000 not in index. {df.shape=}'
            causal_rank = list(df.index).index('causal_00000') + 1
            datapoint['causal_rank'] = causal_rank
            datapoints.append(datapoint)
        except AssertionError as e:
            print(f'{simulation=}: {e}')
            datapoint['causal_rank'] = np.nan
            datapoints.append(datapoint)
        except FileNotFoundError as e:
            print(f'{simulation=}: {e}')
            datapoint['causal_rank'] = np.nan
            datapoints.append(datapoint)

    df = pd.DataFrame(datapoints)

    # rename columns
    df = df.rename(columns={
        'n': 'Number of genomes',
        'e': 'Effect size',
        'r': 'Replicate',
        'causal_rank': 'Rank of causal gene'
    })

    os.makedirs('out', exist_ok=True)
    df.to_csv('out/results.tsv', sep='\t')
    return df


def run_scoary():
    simulations = os.listdir('simulations')
    simulations = list(set([x.split('-')[0] for x in simulations]))
    n_simulations = len(simulations)

    for i, simulation in enumerate(simulations, start=1):
        _scoary(f'{i}/{n_simulations}', simulation)


def plot_all(df: pd.DataFrame, effect_sizes: [float] = [0.5, 1., 1.5, 2, 3.]):
    mpl.use('module://backend_interagg')

    # fill missing values with max_rank + 20
    max_rank = df['Rank of causal gene'].max()
    df = df.fillna(max_rank + 20)

    fig = plt.figure(figsize=(7, 14))

    def add_normal(ax, mean, sd, x, line_color='black', fill_color='red', alpha: float = 0.5):
        # Calculate mean and standard deviation
        y = norm.pdf(x, mean, sd)
        ax.plot(x, y, color=line_color)
        ax.fill_between(x, y, color=fill_color, alpha=alpha)

    for i, effect_size in enumerate(effect_sizes, start=1):
        effect_size_str = str(effect_size).removesuffix('.0')
        # get axes
        ax_title = fig.add_subplot(5, 1, i)
        ax_lineplot = fig.add_subplot(5, 2, i * 2)
        ax_effect_size = fig.add_subplot(5, 2, i * 2 - 1)

        # plot effect size
        _center = effect_size / 2
        x = np.arange(_center - 5, _center + 5, 0.01)
        ax_effect_size.grid(False)
        ax_effect_size.set_xticks([])
        ax_effect_size.set_yticks([])
        ax_effect_size.set_xlabel(f'Distribution of sampled traits')
        add_normal(ax_effect_size, 0, 1, x, fill_color='#a6cee3')
        add_normal(ax_effect_size, effect_size, 1, x, fill_color='#b2df8a')
        # add a dotted line from [0, 0.41] to [effect_size, 0.41]
        ax_effect_size.plot([0, effect_size], [0.41, 0.41], color='black', linestyle='dotted')
        # add a letter d above the dashed line
        ax_effect_size.text(effect_size / 2, 0.435, f'$ d \equal {effect_size_str} \sigma $', horizontalalignment='center', verticalalignment='center')
        ax_effect_size.set_ylim(0, 0.47)

        # set title
        ax_title.set_title(f'Effect size: {effect_size_str}')
        ax_title.grid(False)
        ax_title.axis('off')

        # plot lineplot
        sns.lineplot(
            x='Number of genomes', y='Rank of causal gene', hue='Effect size',
            palette=sns.color_palette(['black'], 1),
            data=df[df['Effect size'] == effect_size],
            ax=ax_lineplot,
            legend=False
        )
        # ax_lineplot.set_ylim(250, 0)
        # make y axis logarithmic
        ax_lineplot.set_yscale('log')
        ax_lineplot.set_yticks([1, 2, 5, 10, 20, 50, 100])
        # ax_lineplot.get_yaxis().tick_right()
        ax_lineplot.get_yaxis().set_label_position("right")
        ax_lineplot.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
        ax_lineplot.set_xticks(df['Number of genomes'].unique())

    plt.tight_layout()
    # plt.show()
    plt.savefig('out/effect_sizes.svg')


if __name__ == '__main__':
    if not os.path.isfile('out/results.tsv'):
        generate_simulations()
        run_scoary()
        df = analyze_scoary_results()
    else:
        df = pd.read_csv('out/results.tsv', sep='\t', index_col=0)

    plot_all(df, effect_sizes=[0.5, 1., 1.5, 2., 3.])

    print('Complete success.')
