# Pangenome Simulator

1) `generate_simulations()`

Script for simulating a pan-genome. Outputs a Roary-like gene_presence_absence.csv and a Traits file.

This script is based on Ola Brynildsrud's [Simulate_pan_genome](https://github.com/AdmiralenOla/Simulate_pan_genome/).

> [!CAUTION]
> Disclaimer: This script is intended for demonstrating the utility of Scoary2 and may or may not be a realistic 
implementation of how bacterial evolution works.

2) `run_scoary()`

Run Scoary2 on the simulated data.

3) `analyze_scoary_results()`

Parse the output of Scoary2 to find the rank o          f the true trait.

Creates. [results.tsv](out%2Fresults.tsv).

If Scoary2 produces no output (no gene left after multiple testing correction) or if the true trait is not in 
the final list of traits, the rank is set to `nan`.

3) `plot_all()`

Plot the results of the analysis.

Creates `out/effect_sizes.png`.

![effect_sizes.svg](out%2Feffect_sizes.svg)