# 3.2.4 — Método de Box (EVOP simplificado) (Python)

Estado: **disponível**.

- Módulo/função: `box_evo.py` — `box_evo(f, x0, Delta, tolx, kmax=10000, verbose=False)` → `BoxResult`
- Exemplo: `ex03_2_4_box.py`, que reproduz os números do deck 3.2.4 (os mesmos números, e a mesma saída, que o exemplo MATLAB).
- Apontamentos: `notes/pt/` (deck 3.2.4).

## Como correr

Só precisa de `numpy`. Na pasta desta secção:

```bash
python ex03_2_4_box.py
```

## O que o exemplo imprime

A tabela das 20 iterações no Himmelblau a partir de (0,0) com Δ = (2,2) e ε_x = 10⁻⁴, os valores de f nos 4 vértices da «iteração à mão», o Rosenbrock a partir de (−1.5, 2) com Δ₀ = (1,1) e o teste das diagonais em (0,0). No fim, a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Iteração à mão: f(0,0) = 170; vértices (−1,−1), (1,−1), (−1,1), (1,1): 170, 146, 130, 106.
- Himmelblau: tabela k = 0..9 igual à do slide; 20 iterações, n_f = 4 × 20 + 1 = 81; termina em (3,2), f = 0.
- Rosenbrock (ε_x = 10⁻⁶, como no caderno da comparação; o deck não indica ε_x): primeira avaliação com f < 10⁻⁴ em n_f = 13 445, com Δ = 6.1·10⁻⁵ (≈ 6·10⁻⁵); no fim, 6474 iterações e n_f = 25 897.
- Diagonais em (0,0) com passo 0.25: f = 4.08 a 11.33 (nenhuma melhora f = 1); f(0.25, 0) = 0.953.
- Última linha: `confere com os slides: sim`.

**Contagens.** `nfev = 1 + 2^n · nit`: f(x₀) conta, e em cada iteração avaliam-se **sempre** os 2ⁿ vértices (convenção do deck). O resultado tem os campos `x`, `fx`, `nit`, `nfev`, `ngev`, `nhev`, `history`, `flag` (e `message`, `cols`, `ftrace`, `Delta`), com as mesmas contagens e a mesma ordem dos vértices que a versão MATLAB (a de `dec2bin`). `ftrace` guarda os valores de f pela ordem de avaliação.
Testado com Python 3.11 e numpy 2.4.
