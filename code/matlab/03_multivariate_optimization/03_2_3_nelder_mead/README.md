# 3.2.3 — Nelder–Mead (MATLAB)

Estado: **disponível**.

- Função: `NelderMead.m` — `[x, fx, info] = NelderMead(f, X0, tolx, tolf, kmax, ftarget, verbose)`. O deck usa `fminsearch`; esta é a implementação da UC, com os coeficientes (α, γ, ρ, σ) = (1, 2, ½, ½) e as regras de aceitação e de encolhimento do pseudocódigo do deck («O algoritmo»).
- Exemplo: `ex03_2_3_nelder_mead.m`, que reproduz os números do deck 3.2.3.
- Apontamentos: `notes/pt/` (deck 3.2.3).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex03_2_3_nelder_mead
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex03_2_3_nelder_mead"`). A ajuda da função: `help NelderMead`.

**Simplex inicial.** `X0` pode ser um ponto x₀ (vetor linha): o simplex é x₀, x₀ + hᵢeᵢ com hᵢ = 0.05·max(1, |x₀ᵢ|), o do deck. Ou uma matriz (n+1) × n com os vértices: é o que usam os exemplos do deck (arestas 1.2 no Himmelblau, 0.5 no Rosenbrock). **Paragem:** f(x_w) − f(x_b) ≤ `tolf` e diâmetro do simplex ≤ `tolx`, medido como maxᵢ‖xᵢ − x_b‖ (os valores por omissão são 10⁻⁶ e 10⁻⁶, os do deck), ou `kmax` iterações (500).

## O que o exemplo imprime

A iteração à mão em x₁² + x₂²; a tabela das iterações do Himmelblau a partir de (0,0) (k = 0..12, com a operação de cada iteração), o resultado final e o outro simplex inicial; o custo no Rosenbrock até f < 10⁻⁴ com os dois simplex iniciais; e a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- À mão: reflexão aceite, x_r = (1, −1), f = 2; novo simplex {(−1,0), (1,−1), (1,2)}, f = 1, 2, 5; 1 avaliação.
- Himmelblau, simplex (0,0), (1.2,0), (0,1.2): tabela k = 0..12 igual à do slide; 54 iterações, n_f = 106, x* = (3, 2), f ≈ 10⁻¹². Com (0,0), (−1.2,0), (0,−1.2): x* = (−2.81; 3.13).
- Rosenbrock de (−1.5, 2): n_f até f < 10⁻⁴ = 166 (arestas 0.5) e 139 (arestas de 5 %).
- Última linha: `confere com os slides: sim`.

**Contagens.** `info.nfev` inclui as n + 1 avaliações do simplex inicial; cada iteração custa 1–2 avaliações, ou n + 2 com encolhimento. `info.nhit` é a primeira avaliação com f < `ftarget` (`Inf` se nunca). É o custo usado na comparação final do capítulo 3. `info.history` tem uma linha por iteração, `[k n_f f(x_b) f(x_w) x_b]`, e `info.ops` a operação de cada linha.
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave.
