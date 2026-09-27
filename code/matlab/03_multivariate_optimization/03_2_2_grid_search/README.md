# 3.2.2 — Pesquisa em grelha (MATLAB)

Estado: **disponível**.

- Função: `GridSearch.m` — `[x, fx, info] = GridSearch(f, a, b, m, verbose)`, com o `meshgrid` do deck lá dentro (n = 2, como no deck).
- Exemplo: `ex03_2_2_grid_search.m`, que reproduz os números do deck 3.2.2.
- Apontamentos: `notes/pt/` (deck 3.2.2).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex03_2_2_grid_search
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex03_2_2_grid_search"`). A ajuda da função: `help GridSearch`.

O exemplo usa também `NelderMead` (3.2.3). Não há cópia da função nesta pasta: o exemplo chama `uc_setup` (ver `code/README.md`), que põe no caminho todas as pastas do código da UC, e usa a função da pasta `03_2_3_nelder_mead`. Corre a partir de qualquer pasta, desde que o repositório esteja clonado inteiro.

## O que o exemplo imprime

Himmelblau em [−5,5]²: os 5 melhores nós da grelha 10 × 10, a grelha fina à volta do melhor nó, e o Nelder–Mead a partir dos 5 melhores nós (ponto inicial, mínimo atingido, n_f), com a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Grelha 10 × 10: 100 avaliações, melhor f = 4.51 em (−2.78, 2.78); os 4 melhores nós, f = 4.51, 4.70, 6.15, 11.88, um em cada bacia.
- Grosseira → fina (mais 100 avaliações em [x − h, x + h], h = 10/9): f = 0.32.
- Nelder–Mead (simplex do deck, tolerâncias 10⁻⁶) a partir dos 4 melhores nós: n_f = 89, 94, 94, 87; os 4 mínimos, f < 10⁻¹¹; total 100 + 364 = 464 avaliações. O 5.º nó leva outra vez a (3, 2), com mais 95.
- Última linha: `confere com os slides: sim`.

Contagens em `info`: `nfev = m²`, `ngev = nhev = 0`, `nit = 0` (uma só passagem pela grelha); `info.history` tem todos os nós ordenados do melhor para o pior (`[x1 x2 f]`) e `info.X1`, `info.X2`, `info.F` a grelha. No Nelder–Mead, `nfev` inclui as 3 avaliações do simplex inicial (o nó da grelha é reavaliado).
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave.
