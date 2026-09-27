# 4.3.2 — Barreira logarítmica (Python)

Estado: **disponível**.

- Módulo/funções: `barrier_log.py` — `barrier_log(f, g, x0, R0, c, tolcomp, tolx, tmax, interno=None, verbose=False)` → `BarrierResult`, com a assinatura e o núcleo do slide «No computador» (P = `pbar(f, g, x, R)`, minimizar a partir de x, u = R/g, parar se J·R ≤ `tolcomp` e ‖x^(t) − x^(t−1)‖/max(1, ‖x^(t−1)‖) ≤ `tolx`, senão R ← cR); `pbar(f, g, x, R, contador=None)`: `inf` se algum g_j ≤ 0 (sem avaliar f), senão f − R Σ ln g_j; com `contador` (dicionário), soma 1 a `contador["fora"]` em cada chamada fora de X.
- Minimizador interno: argumento opcional `interno(P, x)`, que devolve um objeto com `.x` e `.nfev` (chamadas a P) e tem de aceitar P = `inf`. Por omissão, `nelder_mead_comp(P, x)` (tolx = 10⁻¹⁰, tolf = 10⁻¹², kmax = 2000). O slide usa o `fminsearch` do MATLAB.
- `nelder_mead_comp.py` (em `code/python/common/`, o mesmo de 4.3.1; encontrado através de `import uc_setup`) — o Nelder–Mead da UC na versão do script da comparação (ver o README de `04_3_1_exterior_penalty` para as duas diferenças em relação ao `nelder_mead` do deck 3.2.3).
- Exemplo: `ex04_3_2_barrier.py`, que reproduz os números do deck 4.3.2, incluindo a comparação exterior vs. interior (importa `penalty_exterior` da pasta `04_3_1_exterior_penalty/`, no `sys.path` através de `uc_setup`). Os mesmos números que o exemplo MATLAB com `NelderMeadComp`.
- Apontamentos: `notes/pt/` (deck 4.3.2).

## Como correr

Só precisa de `numpy`. Na pasta desta secção:

```bash
python ex04_3_2_barrier.py
```

## O que o exemplo imprime

1. 1D (min x² s.a. x ≥ 1, de x = 2): x_R = 1.366, 1.048, 1.005 para R = 1, 0.1, 0.01.
2. 2D (min x₁² + x₂² s.a. x₁ + x₂ ≥ 1, de (1,1)), R = 1, 0.1, 0.01, 0.001 com arranque a quente: x_R, g(x_R) > 0, ‖x_R − x*‖, ‖x_R − x*‖/R → 1/√2 (erro O(R)), u_R = R/g → 1, e as chamadas a P e fora de X por ciclo.
3. Barreira inversa R/g com R = 0.01: x_R = (0.548; 0.548), contra (0.505; 0.505) da log.
4. Restrição inativa: x_R = (0.882; 0.882) e (0.986; 0.986) para R = 1, 0.1; u_R = R/g → 0.
5. A chamada do slide, com o interno por omissão.
6. A comparação do slide «Exterior vs. interior» e a tabela ciclos / mediana por ciclo / total.

## Saída esperada (resumo)

- Tabelas 1–4 iguais às dos slides.
- Chamada do slide: nit = 10 (R = 10⁻⁹), u = 1.0075, 144–236 chamadas a P por ciclo, `nfev` = 1676 (287 fora de X; n_f efetivo 1389). O slide (com o `fminsearch` do MATLAB) indica nit = 10, u ≈ 1, 140–190 por ciclo. O nit = 10 não depende do minimizador: R no ciclo 9 é 0.1⁸ calculado por multiplicações sucessivas, 1.0000000000000005·10⁻⁸ > `tolcomp` = 10⁻⁸, pelo que J·R ≤ `tolcomp` só pode valer a partir do ciclo 10.
- Comparação (mesmo `nelder_mead_comp`, tolerâncias 10⁻¹⁰): barreira de (1,1) em 5 ciclos, 743 chamadas a P, 26 fora de X, n_f = 717; exterior de (0,0) em 6 ciclos, n_f = 844, e de (1,1), n_f = 864; medianas por ciclo 140 e 146. Com o `nelder_mead` do deck 3.2.3 (contração exterior aceite com ≤): 739 chamadas a P, 26 fora, n_f = 713 (informativo).
- Última linha: `confere com os slides: sim` (o nit da chamada do slide é informativo; o u ≈ 1, a menos de 1 %, entra na verificação).

**Contagens.** Como no slide, `nfev` = chamadas a P + 1 (a avaliação final f(x)): é um **majorante** de n_f, porque as chamadas com x fora de X avaliam g mas não f. O resultado traz também `nPev` (chamadas a P, sem a final), `nfora` (fora de X) e `nfev_efetivo` = nPev − nfora + 1. Os números da comparação (743, 717) são sem a avaliação final, como no slide. O ponto inicial não é avaliado à parte (entra como vértice do simplex do primeiro ciclo). `history` tem uma linha por ciclo t = 0, …, nit: `[t, R, P(ciclo), fora(ciclo), P(acum.), x, u]` (na linha t, R é o valor usado para obter x^(t)); `R` é o do último ciclo; `ngev = nhev = 0`.
Testado com Python 3.11 e numpy 2.4.
