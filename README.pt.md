# Otimização — apontamentos e código

Apontamentos (em português e em inglês) e código MATLAB/Python da unidade curricular **Otimização**
(ISEL — Instituto Superior de Engenharia de Lisboa, Departamento de Matemática).

Autor: **José Firmino Aguilar Madeira** · ORCID [0000-0001-9523-3808](https://orcid.org/0000-0001-9523-3808)

## Conteúdo

| Pasta | O que contém |
|---|---|
| `notes/pt/` | Os 34 decks em PDF, pela ordem da UC: 0 programa, R pré-requisitos, 1–6 capítulos. |
| `notes/en/` | Os mesmos 34 decks em inglês (tradução fiel, para estudantes Erasmus). |
| `references/` | A *Bibliografia da UC* (PDF): `bibliografia_uc.pdf` (português), `course_bibliography.pdf` (inglês). |
| `code/matlab/`, `code/python/` | Uma pasta por método, com o número do deck (p. ex. `02_2_golden_section/`): a função e o exemplo `ex02_2_…` que reproduz os números dos slides. |

## Correr no Colab

Os exemplos em Python correm no Google Colab, sem instalar nada: um caderno por capítulo (ver [`code/README.md`](code/README.md)).

| Capítulo | Colab |
|---|---|
| 2 Otimização unidimensional | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aguilarmadeira/optimization-course/blob/main/code/python/notebooks/02_univariate_optimization.ipynb) |
| 3 Otimização multidimensional sem restrições | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aguilarmadeira/optimization-course/blob/main/code/python/notebooks/03_multivariate_optimization.ipynb) |
| 4 Otimização com restrições | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aguilarmadeira/optimization-course/blob/main/code/python/notebooks/04_constrained_optimization.ipynb) |
| 5 Otimização global | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aguilarmadeira/optimization-course/blob/main/code/python/notebooks/05_global_optimization.ipynb) |
| 6 Otimização multiobjetivo | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aguilarmadeira/optimization-course/blob/main/code/python/notebooks/06_multiobjective_optimization.ipynb) |

## Programa

- **0** Programa e organização da UC · **R** Pré-requisitos e revisão
- **1** Introdução: problemas de otimização e modelação; aplicações
- **2** Otimização unidimensional: secção áurea, interpolação quadrática (Powell), bisseção, Newton–Raphson, derivadas numéricas
- **3** Otimização multidimensional sem restrições: pesquisa aleatória, grelha, Nelder–Mead, Box, Hooke–Jeeves; gradiente, Newton, gradientes conjugados
- **4** Otimização com restrições: Lagrange, KKT, penalização exterior e barreira
- **5** Otimização global: GLODS, *simulated annealing*, algoritmos genéticos
- **6** Otimização multiobjetivo: frente de Pareto, agregação, ε-constrangimento, DMS, NSGA-II

## Código

As funções têm os nomes citados nos slides (p. ex. `GoldenSection.m`); os exemplos começam por `ex` seguido do número do deck, porque MATLAB e Python não aceitam nomes começados por algarismo.
Os exemplos correm a partir de qualquer pasta; a organização do código (utilitários comuns, `uc_setup` para usar as funções nos seus próprios scripts) está em [`code/README.md`](code/README.md).
Os códigos de investigação GLODS e DMS não estão neste repositório: são distribuídos na aula e têm páginas próprias (ver os `README.md` das pastas `05_1_glods` e `06_5_dms`).

## Licenças

- Apontamentos: **CC BY-NC-ND 4.0** (`LICENSE-notes.md`): podem ser partilhados com atribuição, mas não alterados nem usados comercialmente.
- Código: **MIT** (`LICENSE`).

## Como citar

Ver `CITATION.cff` (o GitHub mostra a citação no botão «Cite this repository»).
