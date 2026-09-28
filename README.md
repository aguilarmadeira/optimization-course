# Optimization — lecture notes and code

Lecture notes (in Portuguese and in English) and MATLAB/Python code for the course **Otimização** (Optimization)
at ISEL — Instituto Superior de Engenharia de Lisboa, Department of Mathematics.
Versão em português: [README.pt.md](README.pt.md).

Author: **José Firmino Aguilar Madeira** · ORCID [0000-0001-9523-3808](https://orcid.org/0000-0001-9523-3808)

## Contents

| Folder | Contents |
|---|---|
| `notes/en/` | The 34 slide decks in English (PDF), in course order: 0 syllabus, R prerequisites, chapters 1–6. |
| `notes/pt/` | The same 34 decks in Portuguese (the original version). |
| `references/` | The course bibliography (PDF): `course_bibliography.pdf` (English), `bibliografia_uc.pdf` (Portuguese). |
| `code/matlab/`, `code/python/` | One folder per method, numbered as the deck (e.g. `02_2_golden_section/`): the function and an example `ex02_2_…` that reproduces the numbers on the slides. |

## English version

The English decks are a faithful translation of the Portuguese ones: same slides, same examples and numbers (decimal point instead of decimal comma). The code comments and printed messages are in Portuguese.

## Code

The examples run from any folder; how the code is organised (shared utilities, `uc_setup` to use the functions in your own scripts) is described in [`code/README.md`](code/README.md) (in Portuguese).

## Run in Colab

The Python examples run in Google Colab with nothing to install: one notebook per chapter (comments and printed messages in Portuguese).

| Chapter | Colab |
|---|---|
| 2 Univariate optimization | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aguilarmadeira/optimization-course/blob/main/code/python/notebooks/02_univariate_optimization.ipynb) |
| 3 Unconstrained multivariate optimization | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aguilarmadeira/optimization-course/blob/main/code/python/notebooks/03_multivariate_optimization.ipynb) |
| 4 Constrained optimization | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aguilarmadeira/optimization-course/blob/main/code/python/notebooks/04_constrained_optimization.ipynb) |
| 5 Global optimization | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aguilarmadeira/optimization-course/blob/main/code/python/notebooks/05_global_optimization.ipynb) |
| 6 Multiobjective optimization | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aguilarmadeira/optimization-course/blob/main/code/python/notebooks/06_multiobjective_optimization.ipynb) |

## Syllabus

1. Optimization problems and modelling; applications
2. Univariate optimization: golden section, successive quadratic interpolation (Powell's method), bisection, Newton–Raphson, finite differences
3. Unconstrained multivariate optimization: random search, grid search, Nelder–Mead, Box, Hooke–Jeeves; steepest descent, Newton, conjugate gradients
4. Constrained optimization: Lagrange multipliers, KKT conditions, exterior penalty and barrier methods
5. Global optimization: GLODS, simulated annealing, genetic algorithms
6. Multiobjective optimization: Pareto front, aggregation, ε-constraint, DMS, NSGA-II

The research codes GLODS and DMS are not included; see the `README.md` files in `05_1_glods` and `06_5_dms`.

## Licences

Notes: CC BY-NC-ND 4.0 (`LICENSE-notes.md`). Code: MIT (`LICENSE`).

## Citation

See `CITATION.cff`.
