# 2.3 — Interpolação quadrática sucessiva (Powell) (Python)

Estado: **disponível**.

- Módulo/função: `powell_quadratic.py` — `powell_quadratic(f, x1, Delta, tolx, tolf, kmax=100, verbose=False)` → `PowellResult`
- Exemplo: `ex02_3_powell.py`, que reproduz os números do deck 2.3 (os mesmos números que o exemplo MATLAB).
- Apontamentos: `notes/pt/` (deck 2.3).

## Como correr

Só precisa de `numpy`. Na pasta desta secção:

```bash
python ex02_3_powell.py
```

## O que o exemplo imprime

A tabela das iterações do exemplo-guia (8e¹⁻ˣ + 7 ln x, x₁ = 1, Δ = 0,5, tolerâncias 10⁻³, do 1.º teste de 2019), os erros, as contagens, os dois casos de «Quando a parábola não serve» (x⁴/4 − x²/2) com a salvaguarda por reflexão, e a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Iteração 0: a₁ = −0.619, a₂ = 0.828, x̄ = 1.6237, f(x̄) = 7.6806; 4 iterações (k = 0..3), sem reflexões, n_f = 3 + 4 = 7, x* ≈ 1.6094, f(x*) ≈ 7.6804 (exato 1.6093811). O ponto 1.5 fica no trio em todas as iterações.
- Erros |x̄ − x*|: 1.4e-02, 4.5e-03, 1.9e-04, 4.7e-05.
- Trio (0.3; 0.45; 0.6): a₂ < 0, x̄ = −0.46; trio (0.55; 0.6; 0.65): x̄ = 5.31, f(x̄) = 184.
- Última linha: `confere com os slides: sim`.

O resultado tem os campos `x`, `fx`, `nit`, `nfev`, `ngev`, `nhev`, `history`, `flag` (e `message`), com as mesmas contagens que a versão MATLAB.
Testado com Python 3.11 e numpy 2.4.
