# 4.3.1 — Penalização exterior (Python)

Estado: **disponível**.

- Módulo/função: `penalty_exterior.py` — `penalty_exterior(f, g, h, x0, R0, c, tolviol, tolx, tmax, interno=None, verbose=False)` → `PenaltyResult`, com a assinatura e o núcleo do slide «No computador» (P = f + R(Σ⟨g_j⟩² + Σh_ℓ²), minimizar a partir de x, parar se viol ≤ `tolviol` e ‖x^(t) − x^(t−1)‖/max(1, ‖x^(t−1)‖) ≤ `tolx`, senão R ← cR). Campos `u`, `lam` (o `info.lambda` do MATLAB), `nfev`.
- Minimizador interno: argumento opcional `interno(P, x)`, que devolve um objeto com `.x` e `.nfev` (chamadas a P). Por omissão, o Nelder–Mead da UC na versão da comparação: `nelder_mead_comp(P, x)` (tolx = 10⁻¹⁰, tolf = 10⁻¹², kmax = 2000). O slide usa o `fminsearch` do MATLAB.
- `nelder_mead_comp.py` (em `code/python/common/`, partilhado com 4.3.2; encontrado através de `import uc_setup`) — o Nelder–Mead da UC na versão do script da comparação de 4.3.2 (ver abaixo).
- Exemplo: `ex04_3_1_exterior_penalty.py`, que reproduz os números do deck 4.3.1 e o lado exterior da comparação de 4.3.2 (os mesmos números que o exemplo MATLAB com `NelderMeadComp`).
- Apontamentos: `notes/pt/` (deck 4.3.1).

## Como correr

Só precisa de `numpy`. Na pasta desta secção:

```bash
python ex04_3_1_exterior_penalty.py
```

O exemplo usa também, a título informativo, `nelder_mead` da pasta `03_multivariate_optimization/03_2_3_nelder_mead/` (no `sys.path` através de `uc_setup`).

## Porquê `nelder_mead_comp` e não o `nelder_mead` do deck 3.2.3

Os números da comparação (844; 864 de (1,1)) vêm do Nelder–Mead didático de `make_figs_4comp.py`. Tem os mesmos coeficientes (1, 2, ½, ½), o mesmo simplex inicial (arestas de 5 %) e a mesma paragem que o `nelder_mead` da UC, mas difere em dois pormenores: (1) aceita a contração exterior só se f(x_oc) < f(x_r) (o pseudocódigo do deck 3.2.3 aceita com ≤) — com tolx = 10⁻¹⁰ há empates exatos perto da convergência; (2) calcula x_e = x_c + 2(x_c − x_w) e x_oc = x_c + ½(x_c − x_w) (os mesmos pontos em aritmética exata, arredondamentos diferentes). Com o `nelder_mead` do deck 3.2.3 obtém-se 866 e 873 em vez de 844 e 864 (o exemplo mostra-o). `nelder_mead_comp` é uma cópia adaptada com essas duas diferenças; é igual à da pasta `04_3_2_barrier`.

## O que o exemplo imprime

1. 1D (min x² s.a. x ≥ 1): x_R = ½, 0.909 para R = 1, 10.
2. 2D (min x₁² + x₂² s.a. x₁ + x₂ ≥ 1), R = 1, 10, 100, 1000 com arranque a quente: x_R, g(x_R), ‖x_R − x*‖, R·‖x_R − x*‖ → 1/(2√2) (erro O(1/R)), u_R = −2R⟨g⟩ → 1.
3. Restrição inativa: x_R = (1,1), u_R = 0.
4. Escala das restrições (g̃₂ = 100(1 − x₂)): x₁,R, x₂,R e κ(∇²P) para R = 1, 10, 100.
5. Igualdade (min 2x₁² + x₂² s.a. x₁ + x₂ = 1): x_R, h(x_R), λ_R = −2Rh → 4/3.
6. A chamada do slide, com o interno por omissão.
7. Comparação (lado exterior): tabela dos ciclos de (0,0) e de (1,1) até ‖x^(t) − x*‖ ≤ 10⁻⁴.

## Saída esperada (resumo)

- Tabelas 1–5 iguais às dos slides.
- Chamada do slide: nit = 10 (R = 10⁸), u = 1.0010, 134–257 avaliações por ciclo, `nfev` = 1629. O slide (com o `fminsearch` do MATLAB) indica nit = 10, u ≈ 1, 130–210 por ciclo. Com x_R exatos o critério valeria no ciclo 9; o passo do ciclo 9 é 1.3·10⁻⁶ > 10⁻⁶ por causa do erro do minimizador interno ao longo de g = 0 (P mal condicionada). O nit depende do minimizador interno: é informativo e não entra na verificação.
- Comparação: de (0,0), 6 ciclos, n_f = 844 (mediana 140 por ciclo); de (1,1), 6 ciclos, n_f = 864.
- Última linha: `confere com os slides: sim`.

**Contagens.** Cada chamada a P avalia f uma vez: `nfev` = soma das chamadas a P + 1, a avaliação final f(x). As avaliações de g e h não se contam em n_f. O ponto inicial x⁽⁰⁾ não é avaliado à parte (entra como vértice do simplex do primeiro ciclo). Os n_f da comparação (844, 864) são sem a avaliação final, como no slide: `history[t, 3]`. `history` tem uma linha por ciclo t = 0, …, nit: `[t, R, n_f(ciclo), n_f(acum.), viol, x, u, lambda]` (na linha t, R é o valor usado para obter x^(t)); `R` é o do último ciclo; `ngev = nhev = 0`.
Testado com Python 3.11 e numpy 2.4.
