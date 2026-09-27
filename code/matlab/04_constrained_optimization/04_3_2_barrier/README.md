# 4.3.2 — Barreira logarítmica (MATLAB)

Estado: **disponível**.

- Funções: `BarrierLog.m` — `[x, fx, info] = BarrierLog(f, g, x0, R0, c, tolcomp, tolx, tmax, interno, verbose)`, com a assinatura e o núcleo do slide «No computador» (P = `Pbar(f, g, x, R)`, minimizar a partir de x, u = R./g, parar se J·R ≤ `tolcomp` e ‖x^(t) − x^(t−1)‖/max(1, ‖x^(t−1)‖) ≤ `tolx`, senão R ← cR); `Pbar.m` — `v = Pbar(f, g, x, R)`: `Inf` se algum g_j ≤ 0 (sem avaliar f), senão f − R Σ ln g_j. Devolve `info.nfev`, `info.u`.
- Minimizador interno: argumento opcional `interno`, um handle `[xn, Pn, out] = interno(P, x)` com `out.nfev` = chamadas a P; tem de aceitar P = `Inf`. Por omissão, o do slide: `fminsearch(P, x, opt)` com `opt = optimset('TolX', 1e-10, 'TolFun', 1e-10, 'MaxFunEvals', 2000, 'MaxIter', 2000)`.
- `NelderMeadComp.m` (em `code/matlab/common/`, o mesmo de 4.3.1; posto no caminho por `uc_setup`) — o Nelder–Mead da UC na versão do script da comparação (ver o README de `04_3_1_exterior_penalty` para as duas diferenças em relação ao `NelderMead` do deck 3.2.3). É o interno dos exemplos.
- Exemplo: `ex04_3_2_barrier.m`, que reproduz os números do deck 4.3.2, incluindo a comparação exterior vs. interior (usa `PenaltyExterior` da pasta `04_3_1_exterior_penalty/`, posta no caminho por `uc_setup`).
- Apontamentos: `notes/pt/` (deck 4.3.2).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex04_3_2_barrier
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex04_3_2_barrier"`). A ajuda: `help BarrierLog`, `help Pbar`.

## O que o exemplo imprime

1. 1D (min x² s.a. x ≥ 1, de x = 2): x_R = 1.366, 1.048, 1.005 para R = 1, 0.1, 0.01.
2. 2D (min x₁² + x₂² s.a. x₁ + x₂ ≥ 1, de (1,1)), R = 1, 0.1, 0.01, 0.001 com arranque a quente: x_R, g(x_R) > 0, ‖x_R − x*‖, ‖x_R − x*‖/R → 1/√2 (erro O(R)), u_R = R/g → 1, e as chamadas a P e fora de X por ciclo.
3. Barreira inversa R/g com R = 0.01: x_R = (0.548; 0.548), contra (0.505; 0.505) da log.
4. Restrição inativa: x_R = (0.882; 0.882) e (0.986; 0.986) para R = 1, 0.1; u_R = R/g → 0.
5. A chamada do slide, (a) com o `fminsearch` (informativo) e (b) com `NelderMeadComp`.
6. A comparação do slide «Exterior vs. interior» e a tabela ciclos / mediana por ciclo / total.

## Saída esperada (resumo)

- Tabelas 1–4 iguais às dos slides.
- Chamada do slide com `NelderMeadComp`: nit = 10 (R = 10⁻⁹), u = 1.0075, 144–236 chamadas a P por ciclo, `info.nfev` = 1676 (287 fora de X; n_f efetivo 1389).
- Chamada do slide com o `fminsearch` do **Octave 8.4**: nit = 10 (R = 10⁻⁹), u = 0.9825, 147–176 chamadas a P por ciclo (174 fora de X). O slide indica nit = 10 (R = 10⁻⁹), u ≈ 1, 140–190 por ciclo. O nit = 10 não depende do minimizador: R no ciclo 9 é 0.1⁸ calculado por multiplicações sucessivas, 1.0000000000000005·10⁻⁸ > `tolcomp` = 10⁻⁸, pelo que J·R ≤ `tolcomp` só pode valer a partir do ciclo 10. O u final depende do interno (em R = 10⁻⁹, g(x) ≈ 10⁻⁹ e um erro de 10⁻¹¹ em x muda u em 1–2 %).
- Comparação (mesmo `NelderMeadComp`, tolerâncias 10⁻¹⁰): barreira de (1,1) em 5 ciclos, 743 chamadas a P, 26 fora de X, n_f = 717; exterior de (0,0) em 6 ciclos, n_f = 844, e de (1,1), n_f = 864; medianas por ciclo 140 e 146. Com o `NelderMead` do deck 3.2.3 (contração exterior aceite com ≤): 739 chamadas a P, 26 fora, n_f = 713 (informativo).
- Última linha: `confere com os slides: sim` (o nit da chamada do slide é informativo; o u ≈ 1, a menos de 1 %, entra na verificação com `NelderMeadComp`).

**Contagens.** Como no slide, `info.nfev` = chamadas a P + 1 (a avaliação final f(x)): é um **majorante** de n_f, porque as chamadas com x fora de X avaliam g mas não f. `Pbar` conta essas chamadas numa variável persistente (`Pbar('reset')`, `Pbar('fora')`), e `BarrierLog` devolve também `info.nPev` (chamadas a P, sem a final), `info.nfora` (fora de X) e `info.nfev_efetivo` = nPev − nfora + 1. Os números da comparação (743, 717) são sem a avaliação final, como no slide. O ponto inicial não é avaliado à parte (entra como vértice do simplex do primeiro ciclo). `info.history` tem uma linha por ciclo t = 0, …, nit: `[t R P(ciclo) fora(ciclo) P(acum.) x u]` (na linha t, R é o valor usado para obter x^(t)); `info.R` é o R do último ciclo; `info.ngev = info.nhev = 0`.
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave.
