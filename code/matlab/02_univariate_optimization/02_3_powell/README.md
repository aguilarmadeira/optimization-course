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

A tabela das iterações do exemplo-guia (8e¹⁻ˣ + 7 ln x, x₁ = 1, Δ = 0,5, tolerâncias 10⁻³, do 1.º teste de 2019), os erros, as contagens, os dois casos de «Quando a parábola não serve» (x⁴/4 − x²/2) com a salvaguarda por reflexão, e a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Iteração 0: a₁ = −0.619, a₂ = 0.828, x̄ = 1.6237, f(x̄) = 7.6806; 4 iterações (k = 0..3), sem reflexões, n_f = 3 + 4 = 7, x* ≈ 1.6094, f(x*) ≈ 7.6804 (exato 1.6093811). O ponto 1.5 fica no trio em todas as iterações.
- Erros |x̄ − x*|: 1.4e-02, 4.5e-03, 1.9e-04, 4.7e-05.
- Trio (0.3; 0.45; 0.6): a₂ < 0, x̄ = −0.46; trio (0.55; 0.6; 0.65): x̄ = 5.31, f(x̄) = 184.
- Última linha: `confere com os slides: sim`.

Contagens em `info`: `nfev` (avaliações de f), `ngev` (de f'), `nhev` (de f''), `nit`, e a tabela das iterações em `info.history`.
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave.
