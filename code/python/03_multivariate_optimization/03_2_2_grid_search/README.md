# 3.2.2 — Pesquisa em grelha (Python)

Estado: **disponível**.

- Módulo/função: `grid_search.py` — `grid_search(f, a, b, m, verbose=False)` → `GridSearchResult`, com o `meshgrid` do deck lá dentro (n = 2, como no deck).
- Exemplo: `ex03_2_2_grid_search.py`, que reproduz os números do deck 3.2.2 (os mesmos números que o exemplo MATLAB).
- Apontamentos: `notes/pt/` (deck 3.2.2).

## Como correr

Só precisa de `numpy`. Na pasta desta secção:

```bash
python ex03_2_2_grid_search.py
```

O exemplo usa também `nelder_mead` (3.2.3). Não há cópia do módulo nesta pasta: o exemplo importa `uc_setup` (ver `code/README.md`), que acrescenta ao `sys.path` todas as pastas do código da UC, e usa o módulo da pasta `03_2_3_nelder_mead`. Corre a partir de qualquer pasta, desde que o repositório esteja clonado inteiro.

## O que o exemplo imprime

Himmelblau em [−5,5]²: os 5 melhores nós da grelha 10 × 10, a grelha fina à volta do melhor nó, e o Nelder–Mead a partir dos 5 melhores nós (ponto inicial, mínimo atingido, n_f), com a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Grelha 10 × 10: 100 avaliações, melhor f = 4.51 em (−2.78, 2.78); os 4 melhores nós, f = 4.51, 4.70, 6.15, 11.88, um em cada bacia.
- Grosseira → fina (mais 100 avaliações em [x − h, x + h], h = 10/9): f = 0.32.
- Nelder–Mead (simplex do deck, tolerâncias 10⁻⁶) a partir dos 4 melhores nós: n_f = 89, 94, 94, 87; os 4 mínimos, f < 10⁻¹¹; total 100 + 364 = 464 avaliações. O 5.º nó leva outra vez a (3, 2), com mais 95.
- Última linha: `confere com os slides: sim`.

O resultado tem os campos `x`, `fx`, `nit` (= 0), `nfev` (= m²), `ngev`, `nhev`, `history` (todos os nós, do melhor para o pior: `[x1, x2, f]`), `flag` (e `message`, `h`, `X1`, `X2`, `F`), com as mesmas contagens que a versão MATLAB.
Testado com Python 3.11 e numpy 2.4.
