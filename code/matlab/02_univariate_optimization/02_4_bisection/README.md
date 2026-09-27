# 2.4 — Bisseção (MATLAB)

Estado: **disponível**.

- Função: `Bisection.m` — `[x, fx, info] = Bisection(f, df, a, b, N, verbose)` (`f` só para leitura; pode ser `[]`)
- Exemplo: `ex02_4_bisection.m`, que reproduz os números do deck 2.4.
- Apontamentos: `notes/pt/` (deck 2.4).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex02_4_bisection
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex02_4_bisection"`). A ajuda da função: `help Bisection`.

## O que o exemplo imprime

As tabelas das iterações (k = 0..14) dos exemplos 2 (x² + 54/x em [1,4]) e 1 (−2 sin x + x²/16 em [−1,3]), os valores finais e as contagens, os dois exercícios «Para resolver», o exemplo «Quantas iterações?» e a verificação com os slides. A coluna f(x_k) é só para leitura. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Exemplo 2: x* ≈ 3.0000, b − a = 3/2¹⁵ = 9.2e-05, n_g = N + 2 = 17, n_f = 0; a secção áurea (L = 0.0022) fica 24 vezes maior.
- Exemplo 1: x* ≈ 1.4783, f(x*) ≈ −1.8549, b − a = 1.2e-04.
- Para resolver: 2.5678 / 8.1355 e 1.4496 / −2.7153.
- Última linha: `confere com os slides: sim`.

Contagens em `info`: `nfev` (avaliações de f), `ngev` (de f'), `nhev` (de f''), `nit`, e a tabela das iterações em `info.history`.
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave.
