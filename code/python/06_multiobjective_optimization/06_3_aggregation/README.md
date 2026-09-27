# 6.3 — Soma ponderada e Tchebycheff (Python)

Estado: **disponível**.

- Solver: `scipy.optimize.minimize(method='SLSQP')`, como no slide («custos calculados com SLSQP (SciPy)»). A UC não reimplementa o SLSQP.
- Módulo auxiliar: `aggregation.py`, com duas funções pequenas que fazem o varrimento dos pesos:
  - `weighted_sum(F, X0, bounds, W, cons=(), verbose=False)` → `AggResult`: para cada linha w de `W`, min w·F(x) s.a. `cons` ({'type': 'ineq', 'fun': g}, g ≥ 0) e `bounds`.
  - `weighted_tchebycheff(F, X0, bounds, W, zid, cons=(), verbose=False)` → `AggResult`: min θ s.a. θ − wᵢ(fᵢ(x) − zᵢ^id) ≥ 0 (a reformulação sem o max da nota do slide).
  - `X0`: um ponto → arranque a quente (como no slide); vários → para cada w parte de todos e fica com o melhor.
  - `AggResult`: `x`, `fx`, `nit`, `nfev` (chamadas a F = soma de `r.nfev`, incluindo as dos gradientes por diferenças finitas), `ngev = nhev = 0`, `history` (`[w, x, F(x), n_f]`, e θ no Tchebycheff), `cols`, `status`, `flag`, `message`.
- Exemplo: `ex06_3_aggregation.py`, que reproduz os números do deck 6.3 (os mesmos pontos que o exemplo MATLAB).
- Apontamentos: `notes/pt/` (deck 6.3).

## Como correr

Precisa de `numpy` e `scipy`. Na pasta desta secção:

```bash
python ex06_3_aggregation.py
```

## O que o exemplo faz

1. **Exemplo convexo**, 11 pesos w₁ = 0, 0.1, …, 1, x₀ = (0.5, 0.3), a quente: tabela do slide (x₁* = w₂, f₁ = w₂², f₂ = (1 − w₂)²) e **n_f = 75** (6 + 9 × 7 + 6), como no script das figuras.
2. **Exemplo não convexo**, três pontos iniciais por peso: 11, 100 (solver) e 10 000 pesos (pesquisa exaustiva) → só os 2 extremos com w₁, w₂ > 0. Aviso: com w₁ = 0.5 e um só arranque em (0.5; 0.1), o SLSQP pára em (0.5; 0), um ponto estacionário que é um **máximo** de f̃ em x₁ — por isso os três pontos iniciais.
3. **Tchebycheff ponderado**, z^id = (0, 0), 11 pesos (wᵢ ≥ 10⁻³): pontos no cruzamento w₁f₁ = w₂(1 − f₁²), 9 *non-supported*.
4. **A escala importa**: 81.35 / 85.81 / 100.45 (slide: 81,4 / 85,8 / 100,5), ganha A; normalizado → C / B / A.

## Saída esperada (resumo)

- Convexo: a tabela dos 11 pesos; `n_f = 75`.
- Não convexo: `2 pontos distintos: (0; 1) (1; 0)` com 11, 100 e 10 000 pesos.
- Tchebycheff: f₁ = 1.0000 (exato 0.9995), 0.9460, 0.8828, …, 0.1098, 0.0010.
- Escala: `escolhe C / B / A`.
- Última linha: `confere com os slides: sim`.

As contagens (n_f = 75) são as do SLSQP do SciPy e só se verificam em Python; o MATLAB/Octave usa outro solver e dá os mesmos pontos com outras contagens.
Testado com Python 3.11, numpy 2.4 e SciPy 1.17.
