# 3.3.2 — Newton (puro e amortecido) (MATLAB)

Estado: **disponível**.

- Função: `NewtonND.m` — `[x, fx, info] = NewtonND(f, grad, hess, x0, tolg, kmax, opts)` (`opts.damped`, `opts.modify`, `opts.ls`, `opts.lstol`, `opts.verbose`; no Newton puro `f` pode ser `[]`)
- Pesquisa em linha do amortecido: `LineSearch.m`, em `code/matlab/common/` (a mesma de 3.3.1; ver o README de 3.3.1)
- Exemplo: `ex03_3_2_newton.m`, que reproduz os números do deck 3.3.2.
- Apontamentos: `notes/pt/` (deck 3.3.2).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex03_3_2_newton
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex03_3_2_newton"`). A ajuda da função: `help NewtonND`.

## O que o exemplo imprime

A direção de Newton na quadrática ½(x₁² + 9x₂²) a partir de (9,1); a tabela do Newton puro no Rosenbrock a partir de (−1.5, 2) (k = 0..7) e a do amortecido (k = 0..14); o amortecido com ε_g = 10⁻⁴ (quadro de consulta do 3.3.3); e o Himmelblau a partir de (0,0), puro e amortecido. No fim, a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Quadrática: d₀ = (−9, −1), x₁ = (0, 0): **1 iteração** (n_g = 2, n_H = 1).
- Rosenbrock, puro (ε_g = 10⁻⁸): tabela k = 0..7 igual à do slide; f sobe em k = 2 (751.61) e k = 4 (39.70); **7 iterações**, f = 0 (< 10⁻¹⁶); n_g = 8, n_H = 7.
- Rosenbrock, amortecido: **14 iterações**, f monótona; últimos minimizantes da linha 0.91; 0.98; 1.00. Com ε_g = 10⁻⁴: 14 it., n_g = 15, n_H = 14.
- Himmelblau de (0,0): λ(H₀) = (−42, −26); puro: **4 iterações** para o máximo local (−0.2708, −0.9230), f = 181.6165; amortecido: na 1.ª iteração gᵀd ≥ 0 → usa −g; chega a (3, 2) em **5 iterações**.
- Última linha: `confere com os slides: sim`.

**Puro e amortecido.** Puro (por omissão): resolve H_k d_k = −g_k com `H \ g` e dá o passo completo α = 1. Amortecido (`opts.damped = true`): o núcleo do slide «No computador» — se gᵀd ≥ 0, usa d = −g (salvaguarda) — e α pela pesquisa em linha (`LineSearch`, Brent com tol 10⁻¹⁰, a das figuras). Com `opts.modify = true` faz também o passo 3 do slide «Newton amortecido»: se `chol` falhar, usa H + μI, com μ = 10⁻³·max(1, ‖H‖_F), 10μ, 100μ, … até ser definida positiva (é a variante do caderno da comparação). No Rosenbrock as duas variantes coincidem (H_k é sempre definida positiva no percurso); no Himmelblau de (0,0) a variante com H + μI chega a (3,2) em 4 iterações — as 5 do slide são as da salvaguarda gᵀd ≥ 0 → −g.
Nota: o slide diz «pesquisa em linha … começando em α = 1» (Armijo/Wolfe); a experiência dos slides (e esta função) usa o minimizante da linha, por isso aparecem α > 1 no meio do percurso.

**Contagens.** `info.ngev = nit + 1` (o último é o do teste de paragem), `info.nhev = nit`. `info.nfev` = f em cada iterando (x₀ incluído) + as avaliações das pesquisas em linha (amortecido); `info.nfev_ls` = só as das pesquisas. No Newton puro f só serve para a tabela e para o critério comum da comparação do capítulo; com `f = []`, n_f = 0.
`info.history`: uma linha por k = 0..nit, `[k x_k' f(x_k) ||g_k|| alpha_{k-1} lambda_min(H_{k-1}) dir_{k-1}]` (dir: 0 Newton, 1 H + μI, 2 −g); λ_min é só para leitura.
Arredondamento: `H \ g` usa Cholesky quando H é simétrica definida positiva, o `numpy.linalg.solve` usa LU; a diferença no último bit muda ligeiramente as contagens da pesquisa em linha do amortecido (n_f = 358 aqui, 357 em Python; as iterações e os valores dos slides não mudam).
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave.
