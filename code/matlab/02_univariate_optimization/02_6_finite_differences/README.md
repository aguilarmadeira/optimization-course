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

A conta à mão em x = 2.5 com h = 0.025 (progressiva, regressiva, centrada, f''), a tabela «h cada vez menor?» (h = 0.1 … 10⁻⁸), a tabela do Newton numérico a partir de x₀ = 1, as contagens e a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- x = 2.5, h = 0.025: progressiva −3.5295, centrada −3.64086, f'' ≈ 8.9127; passo de Newton x₁ = 2.9085; com 6 casas f'' daria 8.9136.
- h = 10⁻⁸: f'' calculada = 0.0000 (a verdadeira é 8.912).
- Newton numérico (h = 0.01|x|, ε_g = 10⁻⁴): estagna em 3.00010; para em k = 6 e devolve x₇; n_f = 7 × 3 + 1 = 22. Com h = 10⁻⁴|x| no fim, uma iteração dá 3.0000000.
- Última linha: `confere com os slides: sim`.

Contagens em `info`: `nfev` (avaliações de f), `ngev` (de f'), `nhev` (de f''), `nit`, e a tabela das iterações em `info.history`.
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave.
