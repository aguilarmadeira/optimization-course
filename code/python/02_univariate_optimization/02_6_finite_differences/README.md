# 2.6 — Derivadas numéricas; Newton numérico (Python)

Estado: **disponível**.

- Módulo/função: `num_derivs.py` — `num_derivs(f, x, h=None)` → `(df, ddf, info)`; `newton_numeric.py` — `newton_numeric(f, x0, tolg, kmax, hrel=0.01, verbose=False)` → `NewtonNumericResult`
- Exemplo: `ex02_6_finite_differences.py`, que reproduz os números do deck 2.6 (os mesmos números que o exemplo MATLAB).
- Apontamentos: `notes/pt/` (deck 2.6).

## Como correr

Só precisa de `numpy`. Na pasta desta secção:

```bash
python ex02_6_finite_differences.py
```

## O que o exemplo imprime

A conta à mão em x = 2.5 com h = 0.025 (progressiva, regressiva, centrada, f''), a tabela «h cada vez menor?» (h = 0.1 … 10⁻⁸), a tabela do Newton numérico a partir de x₀ = 1, as contagens e a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- x = 2.5, h = 0.025: progressiva −3.5295, centrada −3.64086, f'' ≈ 8.9127; passo de Newton x₁ = 2.9085; com 6 casas f'' daria 8.9136.
- h = 10⁻⁸: f'' calculada = 0.0000 (a verdadeira é 8.912).
- Newton numérico (h = 0.01|x|, ε_g = 10⁻⁴): estagna em 3.00010; para em k = 6 e devolve x₇; n_f = 7 × 3 + 1 = 22. Com h = 10⁻⁴|x| no fim, uma iteração dá 3.0000000.
- Última linha: `confere com os slides: sim`.

O resultado de `newton_numeric` tem os campos `x`, `fx`, `nit`, `nfev`, `ngev`, `nhev`, `history`, `flag` (e `message`), com as mesmas contagens que a versão MATLAB; `num_derivs` devolve `(df, ddf, info)` com `info['nfev'] = 3`.
Testado com Python 3.11 e numpy 2.4.
