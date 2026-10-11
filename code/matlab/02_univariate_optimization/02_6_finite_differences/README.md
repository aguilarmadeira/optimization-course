# 2.6 — Derivadas numéricas; Newton numérico (MATLAB)

Estado: **disponível**.

- Função: `NumDerivs.m` — `[df, ddf, info] = NumDerivs(f, x, h)`; `NewtonNumeric.m` — `[x, fx, info] = NewtonNumeric(f, x0, tolg, kmax, hrel, verbose)`
- Exemplo: `ex02_6_finite_differences.m`, que reproduz os números do deck 2.6.
- Apontamentos: `notes/pt/` (deck 2.6).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex02_6_finite_differences
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex02_6_finite_differences"`). A ajuda das funções: `help NumDerivs`, `help NewtonNumeric`.

## O que o exemplo imprime

A conta à mão no exemplo-guia, em x = 1.5 com h = 0.015 (progressiva, regressiva, centrada, f''), a tabela «h cada vez menor?» (h = 0.1 … 10⁻⁸), a tabela do Newton numérico a partir de x₀ = 1, as contagens e a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- x = 1.5, h = 0.015: progressiva −0.17255, centrada −0.185605, f'' ≈ 1.7411; passo de Newton x₁ = 1.6066; com 6 casas f'' daria 1.7422, com 4 casas 1.7778.
- h = 10⁻⁸: f'' calculada = 0.0000 (a verdadeira é 1.741).
- Newton numérico (h = 0.01|x|, ε_g = 10⁻⁴): estagna em 1.609407 (x* = 1.609381); para em k = 4 e devolve x₅; n_f = 5 × 3 + 1 = 16. Com h = 10⁻⁴|x| no fim, uma iteração dá 1.6093811.
- Última linha: `confere com os slides: sim`.

Contagens em `info`: `nfev` (avaliações de f), `ngev` (de f'), `nhev` (de f''), `nit`, e a tabela das iterações em `info.history`.
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave.
