# 2.4 — Bisseção (Python)

Estado: **disponível**.

- Módulo/função: `bisection.py` — `bisection(f, df, a, b, N, verbose=False)` → `BisectionResult` (`f` só para leitura; pode ser `None`)
- Exemplo: `ex02_4_bisection.py`, que reproduz os números do deck 2.4 (os mesmos números que o exemplo MATLAB).
- Apontamentos: `notes/pt/` (deck 2.4).

## Como correr

Só precisa de `numpy`. Na pasta desta secção:

```bash
python ex02_4_bisection.py
```

## O que o exemplo imprime

As tabelas das iterações (k = 0..14) dos exemplos 2 (x² + 54/x em [1,4]) e 1 (−2 sin x + x²/16 em [−1,3]), os valores finais e as contagens, os dois exercícios «Para resolver», o exemplo «Quantas iterações?» e a verificação com os slides. A coluna f(x_k) é só para leitura. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Exemplo 2: x* ≈ 3.0000, b − a = 3/2¹⁵ = 9.2e-05, n_g = N + 2 = 17, n_f = 0; a secção áurea (L = 0.0022) fica 24 vezes maior.
- Exemplo 1: x* ≈ 1.4783, f(x*) ≈ −1.8549, b − a = 1.2e-04.
- Para resolver: 2.5678 / 8.1355 e 1.4496 / −2.7153.
- Última linha: `confere com os slides: sim`.

O resultado tem os campos `x`, `fx`, `nit`, `nfev`, `ngev`, `nhev`, `history`, `flag` (e `message`), com as mesmas contagens que a versão MATLAB.
Testado com Python 3.11 e numpy 2.4.
