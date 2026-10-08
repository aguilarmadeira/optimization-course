# 3.3.3 — Gradientes conjugados (Fletcher–Reeves) (MATLAB)

Estado: **disponível**.

- Função: `ConjGrad.m` — `[x, fx, info] = ConjGrad(f, grad, x0, tolg, kmax, opts)` (`opts.restart`, `opts.ls`, `opts.lstol`, `opts.verbose`)
- Pesquisa em linha: `LineSearch.m`, em `code/matlab/common/` (a mesma de 3.3.1; ver o README de 3.3.1)
- Exemplo: `ex03_3_3_conjugate_gradients.m`, que reproduz os números do deck 3.3.3.
- Apontamentos: `notes/pt/` (deck 3.3.3). A tabela comparativa do fim do capítulo está em `../03_3_comparison/`.

**Notação.** Nos slides e nas aulas: S_k é a direção de pesquisa (no código, `d`), λ_k o comprimento do passo (no código, `alpha`) e ∇f_k o gradiente em x_k (no código, `g`).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex03_3_3_conjugate_gradients
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex03_3_3_conjugate_gradients"`). A ajuda da função: `help ConjGrad`.

## O que o exemplo imprime

A tabela das iterações na quadrática ½(x₁² + 9x₂²) a partir de (9,1) (x_k, f, ‖g‖, α, β, d_k), com β₁, d₁, a conjugação d₀ᵀAd₁ = 0 e α₁ = 5/9; o Rosenbrock a partir de (−1.5, 2) com e sem reinício; e, a título informativo, o Rosenbrock com a pesquisa do tipo `fminbnd`. No fim, a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Quadrática: x₁ = (7.2, −0.8); ‖g₁‖² = 103.68, β₁ = 103.68/162 = 0.64; d₁ = (−12.96, 1.44); d₀ᵀAd₁ ≈ 0; α₁ = 5/9; x₂ = (0, 0): **2 iterações, n_g = 3**.
- Rosenbrock (ε_g = 10⁻⁴): FR com reinício a cada n = 2, **36 iterações** (n_g = 37, n_f = 775); sem reinício periódico, **102 iterações**.
- Com a pesquisa do tipo `fminbnd` (`TolX` = 10⁻⁴) e reinício a cada n: não converge em 500 iterações — ilustra o «pode não convergir» do slide «No computador».
- Última linha: `confere com os slides: sim`.

**Reinícios** (os testes 1 e 2 do slide, como no núcleo «No computador»): d = −g quando mod(k+1, n) = 0 (`opts.restart`, por omissão n; 0 desliga) ou quando gᵀd ≥ 0 (sempre ativo). `info.nrestart` conta-os.
**Contagens.** `info.ngev = nit + 1`; `info.nfev` = avaliações das pesquisas em linha (no modo `'brent'` cada pesquisa volta a avaliar φ(0) = f(x_k)) + f em cada iterada (x₀ incluído); `info.nfev_ls` = só as das pesquisas; `info.nhev = 0`.
`info.history`: uma linha por k = 0..nit, `[k x_k' f(x_k) ||g_k|| alpha_{k-1} beta_k d_k']` (β e d da última linha não são calculados: NaN).
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave.
