# 4.3.1 — Penalização exterior (MATLAB)

Estado: **disponível**.

> **Notação das aulas.** Nos decks 4.2 e 4.3 o multiplicador de uma desigualdade g_j ≥ 0 chama-se **λ_j** (e o de uma igualdade **β_ℓ**). No código mantém-se o nome `u` (`info.u`) para as estimativas de λ_j; os valores são os mesmos.

- Função: `PenaltyExterior.m` — `[x, fx, info] = PenaltyExterior(f, g, h, x0, R0, c, tolviol, tolx, tmax, interno, verbose)`, com a assinatura e o núcleo do slide «No computador» (P = f + R(Σ⟨g_j⟩² + Σh_ℓ²), minimizar a partir de x, parar se viol ≤ `tolviol` e ‖x^(t) − x^(t−1)‖/max(1, ‖x^(t−1)‖) ≤ `tolx`, senão fazer R = cR). Devolve `info.nfev`, `info.u`, `info.lambda`.
- Minimizador interno: argumento opcional `interno`, um handle `[xn, Pn, out] = interno(P, x)` com `out.nfev` = chamadas a P. Por omissão, o do slide: `fminsearch(P, x, opt)` com `opt = optimset('TolX', 1e-10, 'TolFun', 1e-10, 'MaxFunEvals', 2000, 'MaxIter', 2000)` e `out.nfev = funcCount`.
- `NelderMeadComp.m` (em `code/matlab/common/`, partilhado com 4.3.2; posto no caminho por `uc_setup`) — o Nelder–Mead da UC na versão do script da comparação de 4.3.2 (ver abaixo). É o interno dos exemplos: `@(P, x) NelderMeadComp(P, x, 1e-10, 1e-12, 2000)`.
- Exemplo: `ex04_3_1_exterior_penalty.m`, que reproduz os números do deck 4.3.1 e o lado exterior da comparação de 4.3.2.
- Apontamentos: `notes/pt/` (deck 4.3.1).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex04_3_1_exterior_penalty
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex04_3_1_exterior_penalty"`). A ajuda: `help PenaltyExterior`, `help NelderMeadComp`. O exemplo usa também, a título informativo, `NelderMead` da pasta `03_multivariate_optimization/03_2_3_nelder_mead/` (no caminho através de `uc_setup`).

## Porquê `NelderMeadComp` e não o `NelderMead` do deck 3.2.3

Os números da comparação (844; 864 de (1,1)) vêm do Nelder–Mead didático de `make_figs_4comp.py`. Tem os mesmos coeficientes (1, 2, ½, ½), o mesmo simplex inicial (arestas de 5 %) e a mesma paragem que o `NelderMead` da UC, mas difere em dois pormenores: (1) aceita a contração exterior só se f(x_oc) < f(x_r) (o pseudocódigo do deck 3.2.3 aceita com ≤) — com tolx = 10⁻¹⁰ há empates exatos perto da convergência; (2) calcula x_e = x_c + 2(x_c − x_w) e x_oc = x_c + ½(x_c − x_w) (os mesmos pontos em aritmética exata, arredondamentos diferentes). Com o `NelderMead` do deck 3.2.3 obtém-se 866 e 873 em vez de 844 e 864 (o exemplo mostra-o). `NelderMeadComp` é uma cópia adaptada com essas duas diferenças; é igual à da pasta `04_3_2_barrier`.

## O que o exemplo imprime

1. 1D (min x² s.a. x ≥ 1): x_R = ½, 0.909 para R = 1, 10.
2. 2D (min x₁² + x₂² s.a. x₁ + x₂ ≥ 1), R = 1, 10, 100, 1000 com arranque a quente: x_R, g(x_R), ‖x_R − x*‖, R·‖x_R − x*‖ → 1/(2√2) (erro O(1/R)), u_R = −2R⟨g⟩ → 1.
3. Restrição inativa: x_R = (1,1), u_R = 0.
4. Escala das restrições (g̃₂ = 100(1 − x₂)): x₁,R, x₂,R e κ(∇²P) para R = 1, 10, 100.
5. Igualdade (min 2x₁² + x₂² s.a. x₁ + x₂ = 1): x_R, h(x_R), β_R = −2Rh → 4/3 (β das aulas; em 4.1, λ* = −4/3).
6. A chamada do slide, (a) com o `fminsearch` (informativo) e (b) com `NelderMeadComp`.
7. Comparação (lado exterior): tabela dos ciclos de (0,0) e de (1,1) até ‖x^(t) − x*‖ ≤ 10⁻⁴.

## Saída esperada (resumo)

- Tabelas 1–5 iguais às dos slides.
- Chamada do slide com `NelderMeadComp`: nit = 10 (R = 10⁸), u = 1.0010, 134–257 avaliações por ciclo, `info.nfev` = 1629.
- Chamada do slide com o `fminsearch` do **Octave 8.4**: nit = 9 (R = 10⁷), u = 1.0000, 156–191 por ciclo. O slide (com o `fminsearch` do MATLAB) indica nit = 10 (R = 10⁸). Não se força: com x_R exatos o critério já vale no ciclo 9 (viol = 5·10⁻⁸, passo 3.2·10⁻⁷ ≤ 10⁻⁶); nit = 10 acontece quando o erro do minimizador interno ao longo de g = 0 (onde P é mal condicionada) faz o passo do ciclo 9 passar 10⁻⁶ — com `NelderMeadComp` é 1.3·10⁻⁶. Depende do minimizador interno e das suas opções, que o slide não mostra.
- Comparação: de (0,0), 6 ciclos, n_f = 844 (mediana 140 por ciclo); de (1,1), 6 ciclos, n_f = 864.
- Última linha: `confere com os slides: sim` (o nit da chamada do slide é informativo e não entra na verificação; o u ≈ 1 entra).

**Contagens.** Cada chamada a P avalia f uma vez: `info.nfev` = soma das chamadas a P (via `out.nfev`/`funcCount`) + 1, a avaliação final f(x). As avaliações de g e h não se contam em n_f. O ponto inicial x⁽⁰⁾ não é avaliado à parte (entra como vértice do simplex do primeiro ciclo). Os n_f da comparação (844, 864) são sem a avaliação final, como no slide: `info.history(t+1, 4)`. `info.history` tem uma linha por ciclo t = 0, …, nit: `[t R n_f(ciclo) n_f(acum.) viol x u lambda]` (na linha t, R é o valor usado para obter x^(t)); `info.R` é o R do último ciclo; `info.ngev = info.nhev = 0`.
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave.
