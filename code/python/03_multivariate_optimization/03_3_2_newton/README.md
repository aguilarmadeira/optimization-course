# 3.3.2 — Newton (puro e amortecido) (Python)

Estado: **disponível**.

- Módulo/função: `newton_nd.py` — `newton_nd(f, grad, hess, x0, tolg=1e-8, kmax=50, damped=False, modify=False, ls="brent", lstol=None, verbose=False)` → `NewtonNDResult` (`ls` = `"brent"` | `"fminbnd"` | `"armijo"`) (no Newton puro `f` pode ser `None`)
- Pesquisa em linha do amortecido: `line_search.py`, em `code/python/common/` (a mesma de 3.3.1; ver o README de 3.3.1)
- Exemplo: `ex03_3_2_newton.py`, que reproduz os números do deck 3.3.2 (os mesmos que o exemplo MATLAB).
- Apontamentos: `notes/pt/` (deck 3.3.2).

**Notação.** Nos slides e nas aulas: S_k é a direção de pesquisa (no código, `d`), λ_k o comprimento do passo (no código, `alpha`) e ∇f_k o gradiente em x_k (no código, `g`). O algoritmo das aulas (S_k = −H_k⁻¹∇f_k e λ_k* ótimo na direção de Newton; numa quadrática com H definida positiva, λ_k* = 1) é `damped=True` com a pesquisa em linha por omissão (`ls="brent"`), que acrescenta a salvaguarda d = −g quando gᵀd ≥ 0; o Newton puro (λ_k = 1) é `damped=False`.

## Como correr

Só precisa de `numpy`. Na pasta desta secção:

```bash
python ex03_3_2_newton.py
```

## O que o exemplo imprime

A direção de Newton na quadrática ½(x₁² + 9x₂²) a partir de (9,1); a tabela do Newton puro no Rosenbrock a partir de (−1.5, 2) (k = 0..7) e a do amortecido (k = 0..14); o amortecido com ε_g = 10⁻⁴ (quadro de consulta do 3.3.3); o amortecido com recuo de Armijo (slide «Que α aceitar?»); e o Himmelblau a partir de (0,0), puro, amortecido e com Armijo (−g ou H + μI). No fim, a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Quadrática: d₀ = (−9, −1), x₁ = (0, 0): **1 iteração** (n_g = 2, n_H = 1).
- Rosenbrock, puro (ε_g = 10⁻⁸): tabela k = 0..7 igual à do slide; f sobe em k = 2 (751.61) e k = 4 (39.70); **7 iterações**, f = 0 (< 10⁻¹⁶); n_g = 8, n_H = 7.
- Rosenbrock, amortecido: **14 iterações**, f monótona; últimos minimizantes da linha 0.91; 0.98; 1.00. Com ε_g = 10⁻⁴: 14 it., n_g = 15, n_H = 14.
- Himmelblau de (0,0): λ(H₀) = (−42, −26); puro: **4 iterações** para o máximo local (−0.2708, −0.9230), f = 181.6165; amortecido: na 1.ª iteração gᵀd ≥ 0 → usa −g; chega a (3, 2) em **5 iterações**.
- Rosenbrock, amortecido com Armijo (α = 1, ½, ¼, …; c₁ = 10⁻⁴): em k = 1, α = 1 dá f = 751.61, ½ dá 49.74, ¼ dá 7.15 (rejeitados) e ⅛ dá 5.24 (aceite); **22 iterações**, n_f = **29**, f monótona, α = 1 nas últimas 7 (n_f ≈ 360 com o minimizante da linha).
- Himmelblau de (0,0) com Armijo: −g na 1.ª iteração, 6 it., n_f = 14; H + μI, 8 it., n_f = 10; ambos em (3, 2).
- Última linha: `confere com os slides: sim`.

**Puro e amortecido.** Puro (por omissão): resolve H_k d_k = −g_k com `np.linalg.solve` e dá o passo completo α = 1. Amortecido (`damped=True`): o núcleo do slide «No computador» — se gᵀd ≥ 0, usa d = −g (salvaguarda) — e α pela pesquisa em linha (`line_search`, Brent com tol 10⁻¹⁰, a das figuras). Com `modify=True` faz também o passo 3 do slide «Newton amortecido»: se a fatorização de Cholesky falhar, usa H + μI, com μ = 10⁻³·max(1, ‖H‖_F), 10μ, 100μ, … até ser definida positiva (é a variante do caderno da comparação). No Rosenbrock as duas variantes coincidem; no Himmelblau de (0,0) a variante com H + μI chega a (3,2) em 4 iterações — as 5 do slide são as da salvaguarda gᵀd ≥ 0 → −g.
**Armijo.** Com `ls="armijo"`, o amortecido não minimiza a linha: tenta α = 1, ½, ¼, … até f(x + αd) ≤ f(x) + c₁α gᵀd (c₁ = 10⁻⁴; no máximo 60 tentativas) e reaproveita o valor aceite como f na nova iterada. Sem a condição de Wolfe. Por omissão (Brent) usa-se o minimizante da linha, por isso aparecem α > 1 no meio do percurso.

**Contagens.** `ngev = nit + 1` (o último é o do teste de paragem), `nhev = nit`. `nfev` = f em cada iterada (x₀ incluído) + as avaliações das pesquisas em linha (amortecido); `nfev_ls` = só as das pesquisas. No Newton puro f só serve para a tabela e para o critério comum da comparação do capítulo; com `f=None`, n_f = 0.
`history`: uma linha por k = 0..nit, `[k, x_k, f(x_k), ||g_k||, alpha_{k-1}, lambda_min(H_{k-1}), dir_{k-1}]` (dir: 0 Newton, 1 H + μI, 2 −g); λ_min é só para leitura.
Arredondamento: `np.linalg.solve` usa LU (como o script das figuras e o caderno), o `H \ g` do MATLAB/Octave usa Cholesky quando H é definida positiva; a diferença no último bit muda ligeiramente as contagens da pesquisa em linha do amortecido (n_f = 357 aqui, 358 em Octave; as iterações e os valores dos slides não mudam).
O resultado tem os campos `x`, `fx`, `nit`, `nfev`, `ngev`, `nhev`, `history`, `flag` (e `message`, `nfev_ls`), com as mesmas contagens que a versão MATLAB.
Testado com Python 3.11 e numpy 2.4.
