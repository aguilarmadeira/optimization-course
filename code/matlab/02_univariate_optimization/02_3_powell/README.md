# 2.3 — Interpolação quadrática sucessiva (Powell) (MATLAB)

Estado: **disponível**.

- Função: `PowellQuadratic.m` — `[x, fx, info] = PowellQuadratic(f, x1, Delta, tolx, tolf, kmax, verbose)`
- Exemplo: `ex02_3_powell.m`, que reproduz os números do deck 2.3.
- Apontamentos: `notes/pt/` (deck 2.3).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex02_3_powell
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex02_3_powell"`). A ajuda da função: `help PowellQuadratic`.

## O que o exemplo imprime

A tabela das iterações do exemplo (x² + 54/x, x₁ = 1, Δ = 1, tolerâncias 10⁻³), os erros, as contagens, os dois casos de «Quando a parábola não serve» (x⁴/4 − x²/2) com a salvaguarda por reflexão, e a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Iteração 0: a₁ = −24, a₂ = 10, x̄ = 2.7, f(x̄) = 27.29; 4 iterações (k = 0..3), sem reflexões, n_f = 3 + 4 = 7, x* = 3, f(x*) = 27.
- Erros |x̄ − 3|: 3.0e-01, 3.8e-02, 1.3e-03, 5.7e-06.
- Trio (0.3; 0.45; 0.6): a₂ < 0, x̄ = −0.46; trio (0.55; 0.6; 0.65): x̄ = 5.31, f(x̄) = 184.
- Última linha: `confere com os slides: sim`.

Contagens em `info`: `nfev` (avaliações de f), `ngev` (de f'), `nhev` (de f''), `nit`, e a tabela das iterações em `info.history`.
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave.
