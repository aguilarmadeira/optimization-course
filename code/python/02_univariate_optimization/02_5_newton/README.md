# 2.5 — Newton–Raphson 1D (Python)

Estado: **disponível**.

- Módulo/função: `newton_1d.py` — `newton_1d(f, df, ddf, x0, tolg, kmax, tolH=1e-12, verbose=False)` → `NewtonResult` (`f` só para leitura; pode ser `None`)
- Exemplo: `ex02_5_newton.py`, que reproduz os números do deck 2.5 (os mesmos números que o exemplo MATLAB).
- Apontamentos: `notes/pt/` (deck 2.5).

## Como correr

Só precisa de `numpy`. Na pasta desta secção:

```bash
python ex02_5_newton.py
```

## O que o exemplo imprime

As tabelas do exemplo-guia (x² + 54/x, x₀ = 2.5), do exemplo de convergência quadrática (x³ − 2x² + 4, x₀ = 1), do exemplo de maximização (8 + 5x − 3x⁴ − 2x⁶, x₀ = 0.5) e da comparação «Os quatro métodos» (x₀ = 1), com as contagens, a classificação do ponto por f'' e a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Exemplo-guia (ε_g = 10⁻⁸; o slide não fixa ε_g): x₁ = 2.9084381; 4 iterações até x* = 3 (erro 2.4e-12); n_g = 5, n_H = 4 (1 + 1 por iteração, mais f'(x₀)), n_f = 0; a 5.ª iteração dá x₅ = 3 exatamente.
- Cúbica: 3 iterações, |f'(x₃)| = 8.1e-04, mínimo em 4/3.
- Maximização: 5 iterações, x* = 0.6617391, f(x*) = 10.565491, f'' < 0 (máximo).
- De x₀ = 1: erro 4.5e-04 após 5 it. (n_g = 6, n_H = 5), 1.8e-15 após 7.
- Última linha: `confere com os slides: sim`.

O resultado tem os campos `x`, `fx`, `nit`, `nfev`, `ngev`, `nhev`, `history`, `flag` (e `message`), com as mesmas contagens que a versão MATLAB.
Testado com Python 3.11 e numpy 2.4.

## Newton salvaguardado (bisseção + Newton)

- Módulo/função: `newton_safeguarded.py` — `newton_safeguarded(f, df, ddf, a, b, x0, tolg, kmax, verbose=False)` → `SafeguardedResult` (`f` só para leitura; pode ser `None`)
- Exemplo: `ex02_5_newton_safeguarded.py`, que reproduz a tabela da página «Newton salvaguardado: bisseção + Newton» do deck 2.5 (os mesmos números que o exemplo MATLAB).

```bash
python ex02_5_newton_safeguarded.py
```

Mantém um intervalo [a, b] com f'(a) < 0 < f'(b). Em cada iteração tenta o passo de Newton x_N = x_k − f'(x_k)/f''(x_k); se x_N sair de (a, b), ou se não der progresso suficiente (|x_N − x_k| maior do que metade do passo de duas iterações antes, o teste de `rtsafe`), faz um passo de bisseção, x_{k+1} = (a + b)/2; depois atualiza o intervalo com o sinal de f'(x_{k+1}). Pára quando |f'(x_{k+1})| < ε_g.

O exemplo mostra primeiro que Newton puro para f'(x) = arctan x a partir de x₀ = 1.5 diverge (1.5 → −1.694 → 2.321 → −5.114) e depois o salvaguardado em [−1, 2] (ε_g = 10⁻¹⁰, a do script das figuras):

| k | x_k | ponto de Newton | x_{k+1} | novo intervalo |
|---|---|---|---|---|
| 0 | 1.5 | −1.694: fora de (a, b) → bisseção | 0.5 | [−1; 0.5] |
| 1 | 0.5 | −0.0796: aceite | −0.0796 | [−0.0796; 0.5] |
| 2 | −0.0796 | 3.4·10⁻⁴: aceite | 3.4·10⁻⁴ | [−0.0796; 3.4·10⁻⁴] |
| 3 | 3.4·10⁻⁴ | −2.5·10⁻¹¹: aceite | −2.5·10⁻¹¹ | [−2.5·10⁻¹¹; 3.4·10⁻⁴] |

Um passo de bisseção e três de Newton, como no slide; 4 iterações, n_g = 7 (f'(a) e f'(b) para verificar o enquadramento, f'(x₀) e uma por iteração), n_H = 4, n_f = 0. Última linha: `confere com os slides: sim`.
`history`: uma linha por iteração, `[k, x_k, x_N, tipo, x_{k+1}, a, b, f'(x_{k+1})]` (tipo 0 = Newton, 1 = bisseção); `nbis` conta as bisseções.
