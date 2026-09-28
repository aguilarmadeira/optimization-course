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

A tabela das iterações do exemplo (x² + 54/x, x₁ = 1, Δ = 1, tolerâncias 10⁻³), os erros, as contagens, os dois casos de «Quando a parábola não serve» (x⁴/4 − x²/2) com a salvaguarda por reflexão, e a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Iteração 0: a₁ = −24, a₂ = 10, x̄ = 2.7, f(x̄) = 27.29; 4 iterações (k = 0..3), sem reflexões, n_f = 3 + 4 = 7, x* = 3, f(x*) = 27.
- Erros |x̄ − 3|: 3.0e-01, 3.8e-02, 1.3e-03, 5.7e-06.
- Trio (0.3; 0.45; 0.6): a₂ < 0, x̄ = −0.46; trio (0.55; 0.6; 0.65): x̄ = 5.31, f(x̄) = 184.
- Última linha: `confere com os slides: sim`.

O resultado tem os campos `x`, `fx`, `nit`, `nfev`, `ngev`, `nhev`, `history`, `flag` (e `message`), com as mesmas contagens que a versão MATLAB.
Testado com Python 3.11 e numpy 2.4.
