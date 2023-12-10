# Benchmark of pair picking

Goal: compare the performance of Scoary vs Scoary2 pair picking algorithms.

## Output

Raw data: [benchmark.tsv](data%2Fbenchmark.tsv)

![benchmark.png](data%2Fbenchmark.png)

**GLM:** ` time ~ n_isolates + n_genes + n_isolates * n_genes`

- `scoary = 0.0006532479935340877 + -8.266041425328844e-07 * n_isolates + -0.00010316416563699979 * n_genes + 2.8076161350353536e-05 * n_isolates * n_genes`
- `scoary2 = 4.729632447019741e-05 + 1.387879503778183e-05 * n_isolates + -2.187177527361452e-06 * n_genes + 6.866970437111624e-07 * n_isolates * n_genes`

Full output: see [benchmark_picking.py](benchmark_picking.py#L308-L352)

![benchmark_with_GLM.png](data%2Fbenchmark_with_GLM.png)

**PySR:** symbolic regression

Operators: `["+", "*", "exp", inv(x)"]`

"Best" scoring models are `constant * n_genes * n_isolates` for Scoary and Scoary2 with the following coefficients:

- Scoary: `2.6693995e-5`
- Scoary2: `8.678912e-7`

Full output: see [benchmark_picking.py](benchmark_picking.py#L396-L426)

![benchmark_with_PySR.png](data%2Fbenchmark_with_PySR.png)

