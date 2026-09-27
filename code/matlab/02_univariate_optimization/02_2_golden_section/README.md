# 2.2 — Secção áurea (MATLAB)

Estado: **disponível**.

- Função: `GoldenSection.m` — `[x, fx, info] = GoldenSection(f, a, b, N, verbose)`
- Exemplo: `ex02_2_golden_section.m`, que reproduz os números do deck 2.2.
- Apontamentos: `notes/pt/` (deck 2.2).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex02_2_golden_section
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex02_2_golden_section"`). A ajuda da função: `help GoldenSection`.

## O que o exemplo imprime

As tabelas das iterações (k = 0..15) dos exemplos 2 (x² + 54/x em [1,4]) e 1 (−3 sin x + x²/8 em [0,4]), os valores finais e as contagens, o exemplo «Quantas reduções?» e a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Exemplo 2: L₁₅ = 0.0022, x* ≈ 3.0002, f(x*) ≈ 27.0000; n_f = N + 2 = 17 para localizar x*, N + 3 = 18 com f no ponto médio final (`info.nfev_loc`, `info.nfev`).
- Exemplo 1: L₁₅ = 0.0029, x* ≈ 1.4489, f(x*) ≈ −2.7153.
- L₀ = 1, ε_x = 0.01: N ≥ 9.57 → N = 10, n_f = 12.
- Última linha: `confere com os slides: sim`.

Contagens em `info`: `nfev` (avaliações de f), `ngev` (de f'), `nhev` (de f''), `nit`, e a tabela das iterações em `info.history`.
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave.
