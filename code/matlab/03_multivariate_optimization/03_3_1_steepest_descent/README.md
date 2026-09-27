# 3.3.1 — Método do gradiente (MATLAB)

Estado: **disponível**.

- Função: `SteepestDescent.m` — `[x, fx, info] = SteepestDescent(f, grad, x0, tolg, kmax, opts)` (`grad = []` usa `GradFD`; `opts.tolx`, `opts.ls`, `opts.lstol`, `opts.verbose`)
- Gradiente numérico: `GradFD.m` — `[g, nfev] = GradFD(f, x, h)` (diferenças centrais, 2n avaliações de f)
- Pesquisa em linha: `LineSearch.m` — `[alpha, nfev, fvals] = LineSearch(phi, metodo, tol, amax)` (em `code/matlab/common/`, partilhada com 3.3.2, 3.3.3 e a comparação; o exemplo põe-na no caminho com `uc_setup`)
- Exemplo: `ex03_3_1_steepest_descent.m`, que reproduz os números do deck 3.3.1 (e o do gradiente numérico do deck 3.3).
- Apontamentos: `notes/pt/` (deck 3.3.1).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex03_3_1_steepest_descent
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex03_3_1_steepest_descent"`). A ajuda: `help SteepestDescent`, `help LineSearch`, `help GradFD`.

## O que o exemplo imprime

A «iteração à mão» na quadrática ½(x₁² + 9x₂²) a partir de (9,1); a tabela do percurso em ziguezague (k = 0..8 e a última linha); o número de iterações e as contagens com as duas pesquisas em linha; a mesma corrida com gradiente numérico; o gradiente de Rosenbrock em (−1.5, 2) por diferenças centrais com h = 10⁻², 10⁻⁴, 10⁻⁶; e o Rosenbrock a partir de (−1.5, 2) com 3000 iterações. No fim, a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Iteração à mão: α₀ = 0.2, x₁ = (7.2, −0.8), f = 28.8, g₁ = (7.2, −7.2), g₁ᵀd₀ ≈ 0.
- Ziguezague (ε_g = 10⁻⁶): tabela k = 0..8 igual à do slide; **74 iterações, n_g = 75**; f_{k+1}/f_k = 0.64 em todas as iterações (a cota de Kantorovich com κ = 9).
- n_f = 1863 = 1788 nas pesquisas em linha (24.2 por pesquisa, «cerca de 24») + 75 iterandos.
- Com a pesquisa do tipo `fminbnd`: as mesmas 74 iterações, **6 avaliações por pesquisa**.
- Gradiente numérico: 74 iterações, n_g = 0, cada gradiente soma 2n = 4 a n_f.
- Rosenbrock em (−1.5, 2): ∇f = (−155, −50); com h = 10⁻⁴ o erro é 6.0·10⁻⁶.
- Rosenbrock (ε_g = 10⁻⁴, k_max = 3000): não converge; x = (0.9844, 0.9690), f = 2.43·10⁻⁴; n_g = 3000, n_f = 60 900.
- Última linha: `confere com os slides: sim`.

## A pesquisa em linha e os «cerca de 6» do `fminbnd`

O slide «No computador» mostra `alpha = fminbnd(phi, 0, amax)` e diz que, no exemplo, o `fminbnd` gasta cerca de 6 avaliações de f por pesquisa, e que a pesquisa de alta precisão das figuras gasta cerca de 24. `LineSearch` tem as duas, sem depender do `fminbnd` (corre igual em MATLAB, Octave e Python):

- `'brent'` (por omissão): enquadramento a partir de (0, 10⁻³) e método de Brent (secção áurea + interpolação parabólica) com tolerância relativa 10⁻¹⁰ — tradução linha a linha da pesquisa usada nos scripts das figuras e na tabela comparativa do capítulo (`scipy.optimize.minimize_scalar`, `method='brent'`, `tol=1e-10`). Dá os mesmos pontos e as mesmas contagens, por isso reproduz as 3000 iterações e n_f = 60 900 do Rosenbrock. Na quadrática: 24.2 avaliações por pesquisa.
- `'fminbnd'`: o algoritmo do `fminbnd` (Brent em [0, a_max], `TolX` = 10⁻⁴), com o alargamento do slide (se α ≈ a_max, repetir com 10·a_max). Na quadrática, φ(α) é uma parábola e cada pesquisa gasta **6** avaliações, como o `fminbnd` do Octave: 3 pontos de secção áurea (0.382, 0.618, 0.236), 1 passo parabólico que acerta em α = 0.2 e 2 pontos a ±tol de α para fechar o intervalo. As 74 iterações não mudam, porque em ambos os casos α ≈ 0.2 com erro muito menor do que o necessário para ‖g‖ < 10⁻⁶.

**Contagens.** `info.ngev = nit + 1` (um gradiente por iteração mais o do teste de paragem; se parar por `kmax`, o gradiente em x_kmax não é calculado, por isso n_g = k_max, como na tabela comparativa). `info.nfev` = avaliações das pesquisas em linha (no modo `'brent'` cada pesquisa volta a avaliar φ(0) = f(x_k), a primeira avaliação do enquadramento) + f em cada iterando (x₀ incluído); `info.nfev_ls` = só as das pesquisas em linha. Com `grad = []`, cada gradiente custa 2n avaliações de f (em `nfev`) e `ngev = 0`. `info.nhev = 0`.
`info.history`: uma linha por k = 0..nit, `[k x_k' f(x_k) ||g_k|| alpha_{k-1}]`.
Paragem: ‖g_k‖ ≤ ε_g (flag 0), k ≥ k_max (flag 1) ou ‖x_{k+1} − x_k‖/max(1, ‖x_k‖) ≤ ε_x (flag 2; `opts.tolx`, por omissão 0).
Nos exemplos a quadrática escreve-se ½(x₁² + 9x₂²) em vez de ½xᵀAx: o produto matricial arredonda de maneira diferente em MATLAB e em numpy e muda o número de avaliações de Brent (24.8 em vez de 24.2 por pesquisa); com a forma escalar as duas versões dão os mesmos números.
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave.
