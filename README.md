# Optimization — lecture notes and code

Lecture notes (in Portuguese) and MATLAB/Python code for the course **Otimização** (Optimization)
at ISEL — Instituto Superior de Engenharia de Lisboa, Department of Mathematics.
Versão em português: [README.pt.md](README.pt.md).

Author: **José Firmino Aguilar Madeira** · ORCID [0000-0001-9523-3808](https://orcid.org/0000-0001-9523-3808)

## Contents

| Folder | Contents |
|---|---|
| `notes/pt/` | The 34 slide decks (PDF, Portuguese), in course order: 0 syllabus, R prerequisites, chapters 1–6. |
| `references/` | The course bibliography (PDF). |
| `code/matlab/`, `code/python/` | One folder per method, numbered as the deck (e.g. `02_2_golden_section/`): the function and an example `ex02_2_…` that reproduces the numbers on the slides. |

## Syllabus

1. Optimization problems and modelling; applications
2. Univariate optimization: golden section, successive quadratic interpolation (Powell), bisection, Newton–Raphson, finite differences
3. Unconstrained multivariate optimization: random search, grid search, Nelder–Mead, Box, Hooke–Jeeves; steepest descent, Newton, conjugate gradients
4. Constrained optimization: Lagrange multipliers, KKT conditions, exterior penalty and barrier methods
5. Global optimization: GLODS, simulated annealing, genetic algorithms
6. Multiobjective optimization: Pareto front, aggregation, ε-constraint, DMS, NSGA-II

The research codes GLODS and DMS are not included; see the `README.md` files in `05_1_glods` and `06_5_dms`.

## Licences

Notes: CC BY-NC-ND 4.0 (`LICENSE-notes.md`). Code: MIT (`LICENSE`).

## Citation

See `CITATION.cff`.
