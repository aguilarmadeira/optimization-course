# 6.3 — Soma ponderada e Tchebycheff (MATLAB)

Estado: **disponível**.

- O deck não tem uma função própria da UC: o slide do 6.4 mostra a soma ponderada como `fmincon(@(x) w1*f1(x)+w2*f2(x), ...)` sem `nonl`. Duas funções auxiliares pequenas fazem o varrimento dos pesos:
  - `WeightedSum.m` — `[x, fx, info] = WeightedSum(F, X0, lb, ub, g, W, verbose)`: para cada linha w de `W`, min w·F(x) s.a. g(x) ≥ 0 (forma da UC), lb ≤ x ≤ ub.
  - `WeightedTchebycheff.m` — `[x, fx, info] = WeightedTchebycheff(F, X0, lb, ub, g, W, zid, verbose)`: min maxᵢ wᵢ(fᵢ − zᵢ^id) pela reformulação da nota do slide, min θ s.a. θ − wᵢ(fᵢ(x) − zᵢ^id) ≥ 0, com z = [x θ].
  - `F` devolve o vetor linha [f₁(x) … f_m(x)]; `X0` tem um ponto por linha. Com um só ponto, arranque a quente (como no slide); com vários, para cada w parte de todos e fica com o melhor (problemas não convexos).
  - `info`: `nfev` (chamadas a F, incluindo as dos gradientes por diferenças finitas do solver; com `fmincon` é a soma de `out.funcCount`), `ngev = nhev = 0`, `nit`, `history` (uma linha por w: `[w x F(x) n_f]`, e θ no Tchebycheff), `cols`, `exitflag`, `solver`, `flag`, `message`. `help WeightedSum`.
- Exemplo: `ex06_3_aggregation.m`, que reproduz os números do deck 6.3.
- Solver: `fmincon` (Optimization Toolbox), se existir; em GNU Octave, `sqp` (Octave base). Os pontos coincidem. Com um objetivo côncavo, o `sqp` do Octave pode parar com erro («QP subproblem is non-convex and unbounded»); as funções apanham o erro e contam esse arranque como falhado (os outros pontos iniciais resolvem).
- Apontamentos: `notes/pt/` (deck 6.3).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex06_3_aggregation
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex06_3_aggregation"`; demora uns 8 s).

## O que o exemplo faz

1. **Exemplo convexo** (f₁ = x₁² + x₂², f₂ = (x₁ − 1)² + x₂², x ∈ [0,1]², x₁² + x₂² ≤ 1): 11 pesos w₁ = 0, 0.1, …, 1, x₀ = (0.5, 0.3), a quente. Tabela do slide: x₁* = w₂, f₁ = w₂², f₂ = (1 − w₂)²; os três pesos de «O que fazem os pesos?».
2. **Exemplo não convexo** (f₁ = x₁, f₂ = 1 − x₁² + x₂, x ∈ [0,1] × [0; 0.6]), três pontos iniciais por peso: w₁ = 0.2, 0.5, 0.8 dão sempre um extremo; 11 e 100 pesos (solver) e 10 000 pesos (pesquisa exaustiva em x₁) → só 2 pontos com w₁, w₂ > 0.
3. **Tchebycheff ponderado** no exemplo não convexo, z^id = (0, 0), 11 pesos (wᵢ ≥ 10⁻³, como no script das figuras): os pontos ficam no cruzamento w₁f₁ = w₂(1 − f₁²), 9 deles interiores (*non-supported*).
4. **A escala importa**: A, B, C; w = (0.9; 0.1) sem normalizar → 81.35 / 85.81 / 100.45 (slide: 81,4 / 85,8 / 100,5), ganha A; normalizado com z^id = (0.5; 800), ẑ^nad = (1.5; 1000) → C / B / A.

## Contagens

O slide dá **n_f = 75** para o exemplo convexo (11 pontos, SLSQP do SciPy, gradientes por diferenças finitas, arranque a quente). Esse número é do SLSQP e **só se verifica na versão Python**. Aqui o solver é outro (`sqp` no Octave: n_f = 106 chamadas a F; `fmincon`: a soma de `out.funcCount`), por isso as contagens imprimem-se só a título informativo e não entram na verificação. Os pontos da frente coincidem. O ponto inicial não é avaliado à parte (só pelo solver); as avaliações finais F(x) para a tabela não contam.

## Saída esperada (resumo)

- Convexo: a tabela dos 11 pesos do slide (x₁* = 1.000, 0.900, …, 0.000; f₁ = 1.000, 0.810, …; f₂ = 0.000, 0.010, …); Octave `sqp`: n_f = 106 (informativo).
- Não convexo: `11 pesos … 2 pontos distintos: (0; 1) (1; 0)`, e o mesmo com 100 e 10 000 pesos.
- Tchebycheff: f₁ = 0.9995, 0.9460, 0.8828, 0.8084, 0.7208, 0.6180, 0.5000, 0.3699, 0.2361, 0.1098, 0.0010.
- Escala: `escolhe C / B / A`.
- Última linha: `confere com os slides: sim`.

Testado em GNU Octave 8.4 (ramo `sqp`); o ramo `fmincon` usa só a sintaxe documentada mas não foi corrido aqui (sem MATLAB). Só usa funções comuns a MATLAB e Octave.
