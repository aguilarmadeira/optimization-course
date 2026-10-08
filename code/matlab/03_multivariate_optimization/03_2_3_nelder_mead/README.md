# 3.2.3 — Nelder–Mead (MATLAB)

Estado: **disponível**.

- Função: `NelderMead.m` — `[x, fx, info] = NelderMead(f, X0, tolx, tolf, kmax, ftarget, verbose, opts)`, com `opts.gamma`, `opts.beta`, `opts.eps` e `opts.regra`. O deck usa `fminsearch`; esta é a implementação da UC.
- Exemplo: `ex03_2_3_nelder_mead.m`, que reproduz os números do deck 3.2.3.
- Apontamentos: `notes/pt/` (deck 3.2.3).

O algoritmo é o das aulas (Deb, 2012, sec. 3.3.2), com a notação dos slides: xl é o melhor vértice, xg o seguinte ao pior, xh o pior e xc o centroide dos vértices exceto xh. Reflexão xr = 2xc − xh; expansão (1+γ)xc − γxh; contração exterior (1+β)xc − βxh; contração interior (1−β)xc + βxh; por omissão γ = 2 e β = 0.5. Nos casos-limite seguem-se as regras do Nelder–Mead padrão (`opts.regra = 'padrao'`): a expansão só fica se f(xe) < f(xr); uma contração só fica se melhorar, senão encolhe-se o simplex para metade em direção a xl. Com `opts.regra = 'deb'` o ponto novo substitui sempre xh e não há encolhimento.

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex03_2_3_nelder_mead
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex03_2_3_nelder_mead"`). A ajuda da função: `help NelderMead`.

**Simplex inicial.** `X0` pode ser um ponto x₀ (vetor linha): o simplex é x₀, x₀ + hᵢeᵢ com hᵢ = 0.05·max(1, |x₀ᵢ|). Ou uma matriz (n+1) × n com os vértices: é o que usam os exemplos do deck (arestas 1.2 no Himmelblau, 0.5 no Rosenbrock).

**Paragem.** Com `opts.eps` (a das aulas): análise do erro Q = √(Σᵢ (f(xᵢ) − f(xc))² / (n+1)) ≤ ε, calculada em cada iteração com o novo simplex e o centroide dessa iteração; f(xc) conta como uma avaliação. Sem `opts.eps` (Nelder–Mead padrão, como o `fminsearch`): f(xh) − f(xl) ≤ `tolf` e maxᵢ‖xᵢ − xl‖ ≤ `tolx` (10⁻⁶ e 10⁻⁶ por omissão). Em ambos os casos, no máximo `kmax` iterações (500).

## O que o exemplo imprime

A iteração à mão em x₁² + x₂²; a tabela das iterações do Himmelblau a partir de (0,0) com ε = 10⁻⁶ (k = 0..12, com Q e a operação de cada iteração), o resultado final e o outro simplex inicial; o custo no Rosenbrock até f < 10⁻⁴ com os dois simplex iniciais, com e sem as avaliações de f(xc); e a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- À mão: reflexão aceite, xr = (1, −1), f = 2; novo simplex {(−1,0), (1,−1), (1,2)}, f = 1, 2, 5; Q = 2.380; 2 avaliações (xr e xc).
- Himmelblau, simplex (0,0), (1.2,0), (0,1.2), ε = 10⁻⁶: tabela k = 0..12 igual à do slide (k = 1: expansão, Q = 57.286; k = 2: reflexão, Q = 52.768); 31 iterações, n_f = 92, x* = (3, 2). Com (0,0), (−1.2,0), (0,−1.2): x* = (−2.81; 3.13).
- Rosenbrock de (−1.5, 2): n_f até f < 10⁻⁴ = 244 (arestas 0.5) e 211 (arestas de 5 %); sem as avaliações de f(xc) (Nelder–Mead padrão): 166 e 139.
- Última linha: `confere com os slides: sim`.

**Contagens.** `info.nfev` inclui as n + 1 avaliações do simplex inicial; cada iteração custa 1–2 avaliações (mais 1 para f(xc) com `opts.eps`), ou mais n com encolhimento. `info.nhit` é a primeira avaliação com f < `ftarget` (`Inf` se nunca). É o custo usado na comparação final do capítulo 3. `info.history` tem uma linha por k = 0, …, nit, `[k n_f f(xl) f(xh) Q xl]` (Q = NaN sem `opts.eps`), e `info.ops` a operação de cada linha.
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave.
