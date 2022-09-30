![Scoary 2 logo (light mode)](media/scoary-2-logo-full.svg#gh-light-mode-only)
![Scoary 2 logo (dark mode)](media/scoary-2-logo-full-dark.svg#gh-dark-mode-only)

Scoary 2 associates orthogenes (e.g. generated using [OrthoFinder][orthofinder]
or [Roary][roary] to traits. It reports a list of genes sorted by strength of
association per trait. The results can be explored interactively with a simple, static HTML/JS app.

# Wiki

- [Home](https://github.com/MrTomRod/scoary-2/wiki/Home)
- [Installation](https://github.com/MrTomRod/scoary-2/wiki/Installation)
- [Usage](https://github.com/MrTomRod/scoary-2/wiki/Usage)
- [Input](https://github.com/MrTomRod/scoary-2/wiki/Input)
- [Output](https://github.com/MrTomRod/scoary-2/wiki/Output)
- [Tutorial](https://github.com/MrTomRod/scoary-2/wiki/Tutorial)
- [App](https://github.com/MrTomRod/scoary-2/wiki/App)
- [Understanding the p values](https://github.com/MrTomRod/scoary-2/wiki/Understanding-the-p-values)
- [Usage as Python library](https://github.com/MrTomRod/scoary-2/wiki/Usage-as-Python-library)

# Todo:

- [X] Binarize traits during multiprocessing
- [X] Idea: min_qval * min_empirical_p as new score
- [X] Improve messaging and logging
- [X] Remove multiple_testing_picking, more like Scoary 1 syntax?
- [X] Log runtime?
- [X] Refactoring
- [X] `overview.html`: go to function, popover always left
- [X] `overview.html` & `trait.html`: ~~set proper title~~, add logo/favicon
- [ ] `overview.html` & `trait.html`: link to `summary.tsv` / `coverage-matrix.tsv` / `result.tsv` download
- [X] Create wiki


[orthofinder]: https://github.com/davidemms/OrthoFinder/
[roary]: https://sanger-pathogens.github.io/Roary/
