# 2.2 — Secção áurea (Python)

Estado: **disponível**.

- Módulo/função: `golden_section.py` — `golden_section(f, a, b, N, verbose=False)` → `GoldenResult`
- Exemplo: `ex02_2_golden_section.py`, que reproduz os números do deck 2.2 (os mesmos números que o exemplo MATLAB).
- Apontamentos: `notes/pt/` (deck 2.2).

## Como correr

Só precisa de `numpy`. Na pasta desta secção:

```bash
python ex02_2_golden_section.py
```

## O que o exemplo imprime

As tabelas das iterações (k = 0..15) dos exemplos 2 (x² + 54/x em [1,4]) e 1 (−3 sin x + x²/8 em [0,4]), os valores finais e as contagens, o exemplo «Quantas reduções?» e a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Exemplo 2: L₁₅ = 0.0022, x* ≈ 3.0002, f(x*) ≈ 27.0000; n_f = N + 2 = 17 para localizar x*, N + 3 = 18 com f no ponto médio final (`nfev_loc`, `nfev`).
- Exemplo 1: L₁₅ = 0.0029, x* ≈ 1.4489, f(x*) ≈ −2.7153.
- L₀ = 1, ε_x = 0.01: N ≥ 9.57 → N = 10, n_f = 12.
- Última linha: `confere com os slides: sim`.

O resultado tem os campos `x`, `fx`, `nit`, `nfev`, `ngev`, `nhev`, `history`, `flag` (e `message`), com as mesmas contagens que a versão MATLAB.
Testado com Python 3.11 e numpy 2.4.
