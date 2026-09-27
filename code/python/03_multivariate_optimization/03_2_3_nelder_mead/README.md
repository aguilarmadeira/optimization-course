# 3.2.3 — Nelder–Mead (Python)

Estado: **disponível**.

- Módulo/função: `nelder_mead.py` — `nelder_mead(f, X0, tolx=1e-6, tolf=1e-6, kmax=500, ftarget=-inf, verbose=False)` → `NelderMeadResult`. O deck mostra `scipy.optimize.minimize(f, x0, method='Nelder-Mead')`; esta é a implementação da UC, com os coeficientes (α, γ, ρ, σ) = (1, 2, ½, ½) e as regras de aceitação e de encolhimento do pseudocódigo do deck («O algoritmo»). Só usa `numpy`.
- Exemplo: `ex03_2_3_nelder_mead.py`, que reproduz os números do deck 3.2.3 (os mesmos números que o exemplo MATLAB).
- Apontamentos: `notes/pt/` (deck 3.2.3).

## Como correr

Só precisa de `numpy`. Na pasta desta secção:

```bash
python ex03_2_3_nelder_mead.py
```

**Simplex inicial.** `X0` pode ser um ponto x₀: o simplex é x₀, x₀ + hᵢeᵢ com hᵢ = 0.05·max(1, |x₀ᵢ|), o do deck. Ou uma matriz (n+1) × n com os vértices: é o que usam os exemplos do deck (arestas 1.2 no Himmelblau, 0.5 no Rosenbrock). **Paragem:** f(x_w) − f(x_b) ≤ `tolf` e diâmetro do simplex ≤ `tolx`, medido como maxᵢ‖xᵢ − x_b‖, ou `kmax` iterações.

## O que o exemplo imprime

A iteração à mão em x₁² + x₂²; a tabela das iterações do Himmelblau a partir de (0,0) (k = 0..12, com a operação de cada iteração), o resultado final e o outro simplex inicial; o custo no Rosenbrock até f < 10⁻⁴ com os dois simplex iniciais; e a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- À mão: reflexão aceite, x_r = (1, −1), f = 2; novo simplex {(−1,0), (1,−1), (1,2)}, f = 1, 2, 5; 1 avaliação.
- Himmelblau, simplex (0,0), (1.2,0), (0,1.2): tabela k = 0..12 igual à do slide; 54 iterações, n_f = 106, x* = (3, 2), f ≈ 10⁻¹². Com (0,0), (−1.2,0), (0,−1.2): x* = (−2.81; 3.13).
- Rosenbrock de (−1.5, 2): n_f até f < 10⁻⁴ = 166 (arestas 0.5) e 139 (arestas de 5 %).
- Última linha: `confere com os slides: sim`.

O resultado tem os campos `x`, `fx`, `nit`, `nfev`, `nhit`, `ngev`, `nhev`, `history`, `ops`, `simplex`, `flag` (e `message`), com as mesmas contagens que a versão MATLAB: `nfev` inclui as n + 1 avaliações do simplex inicial; 1–2 por iteração, n + 2 com encolhimento.
Nota: o `nelder_mead` do caderno `cap3_comparacao.ipynb` é uma versão compacta com outra regra de aceitação das contrações; no Rosenbrock (arestas 0.5) dá o mesmo n_f = 166 até f < 10⁻⁴.
Testado com Python 3.11 e numpy 2.4.
