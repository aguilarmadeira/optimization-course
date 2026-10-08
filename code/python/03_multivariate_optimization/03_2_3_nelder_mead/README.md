# 3.2.3 — Nelder–Mead (Python)

Estado: **disponível**.

- Módulo/função: `nelder_mead.py` — `nelder_mead(f, X0, tolx=1e-6, tolf=1e-6, kmax=500, ftarget=-inf, verbose=False, gamma=2.0, beta=0.5, eps=None, regra="padrao")` → `NelderMeadResult`. O deck mostra `scipy.optimize.minimize(f, x0, method='Nelder-Mead')`; esta é a implementação da UC. Só usa `numpy`.
- Exemplo: `ex03_2_3_nelder_mead.py`, que reproduz os números do deck 3.2.3 (os mesmos números que o exemplo MATLAB).
- Apontamentos: `notes/pt/` (deck 3.2.3).

O algoritmo é o das aulas (Deb, 2012, sec. 3.3.2), com a notação dos slides: x_l é o melhor vértice, x_g o seguinte ao pior, x_h o pior e x_c o centroide dos vértices exceto x_h. Reflexão x_r = 2x_c − x_h; expansão (1+γ)x_c − γx_h; contração exterior (1+β)x_c − βx_h; contração interior (1−β)x_c + βx_h; por omissão γ = 2 e β = 0.5. Nos casos-limite seguem-se as regras do Nelder–Mead padrão (`regra = "padrao"`): a expansão só fica se f(x_e) < f(x_r); uma contração só fica se melhorar, senão encolhe-se o simplex para metade em direção a x_l. Com `regra = "deb"` o ponto novo substitui sempre x_h e não há encolhimento.

## Como correr

Só precisa de `numpy`. Na pasta desta secção:

```bash
python ex03_2_3_nelder_mead.py
```

**Simplex inicial.** `X0` pode ser um ponto x₀: o simplex é x₀, x₀ + hᵢeᵢ com hᵢ = 0.05·max(1, |x₀ᵢ|). Ou uma matriz (n+1) × n com os vértices: é o que usam os exemplos do deck (arestas 1.2 no Himmelblau, 0.5 no Rosenbrock).

**Paragem.** Com `eps` (a das aulas): análise do erro Q = √(Σᵢ (f(xᵢ) − f(x_c))² / (n+1)) ≤ ε, calculada em cada iteração com o novo simplex e o centroide dessa iteração; f(x_c) conta como uma avaliação. Sem `eps` (Nelder–Mead padrão, como o `fminsearch`): f(x_h) − f(x_l) ≤ `tolf` e maxᵢ‖xᵢ − x_l‖ ≤ `tolx`. Em ambos os casos, no máximo `kmax` iterações.

## O que o exemplo imprime

A iteração à mão em x₁² + x₂²; a tabela das iterações do Himmelblau a partir de (0,0) com ε = 10⁻⁶ (k = 0..12, com Q e a operação de cada iteração), o resultado final e o outro simplex inicial; o custo no Rosenbrock até f < 10⁻⁴ com os dois simplex iniciais, com e sem as avaliações de f(x_c); e a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- À mão: reflexão aceite, x_r = (1, −1), f = 2; novo simplex {(−1,0), (1,−1), (1,2)}, f = 1, 2, 5; Q = 2.380; 2 avaliações (x_r e x_c).
- Himmelblau, simplex (0,0), (1.2,0), (0,1.2), ε = 10⁻⁶: tabela k = 0..12 igual à do slide (k = 1: expansão, Q = 57.286; k = 2: reflexão, Q = 52.768); 31 iterações, n_f = 92, x* = (3, 2). Com (0,0), (−1.2,0), (0,−1.2): x* = (−2.81; 3.13).
- Rosenbrock de (−1.5, 2): n_f até f < 10⁻⁴ = 244 (arestas 0.5) e 211 (arestas de 5 %); sem as avaliações de f(x_c) (Nelder–Mead padrão): 166 e 139.
- Última linha: `confere com os slides: sim`.

O resultado tem os campos `x`, `fx`, `nit`, `nfev`, `nhit`, `ngev`, `nhev`, `history`, `ops`, `simplex`, `Q`, `flag` (e `message`), com as mesmas contagens que a versão MATLAB: `nfev` inclui as n + 1 avaliações do simplex inicial; 1–2 por iteração (mais 1 para f(x_c) com `eps`), n mais com encolhimento. `history` tem uma linha por k = 0, …, nit: [k, n_f, f(x_l), f(x_h), Q, x_l].
Nota: o `nelder_mead` do caderno `cap3_comparacao.ipynb` é uma versão compacta da versão das aulas; no Rosenbrock (arestas 0.5) dá o mesmo n_f = 244 até f < 10⁻⁴.
Testado com Python 3.11 e numpy 2.4.
