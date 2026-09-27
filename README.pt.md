# Otimização — apontamentos e código

Apontamentos (em português) e código MATLAB/Python da unidade curricular **Otimização**
(ISEL — Instituto Superior de Engenharia de Lisboa, Departamento de Matemática).

Autor: **José Firmino Aguilar Madeira** · ORCID [0000-0001-9523-3808](https://orcid.org/0000-0001-9523-3808)

## Conteúdo

| Pasta | O que contém |
|---|---|
| `notes/pt/` | Os 34 decks em PDF, pela ordem da UC: 0 programa, R pré-requisitos, 1–6 capítulos. |
| `references/` | A *Bibliografia da UC* (PDF). |
| `code/matlab/`, `code/python/` | Uma pasta por método, com o número do deck (p. ex. `02_2_golden_section/`): a função e o exemplo `ex02_2_…` que reproduz os números dos slides. |

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
Os códigos de investigação GLODS e DMS não estão neste repositório: são distribuídos na aula e têm páginas próprias (ver os `README.md` das pastas `05_1_glods` e `06_5_dms`).

## Licenças

- Apontamentos: **CC BY-NC-ND 4.0** (`LICENSE-notes.md`): podem ser partilhados com atribuição, mas não alterados nem usados comercialmente.
- Código: **MIT** (`LICENSE`).

## Como citar

Ver `CITATION.cff` (o GitHub mostra a citação no botão «Cite this repository»).
