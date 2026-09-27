# 6.6 — NSGA-II (MATLAB)

Estado: **disponível**.

- Função: `NSGA2.m` — `[X, F, info] = NSGA2(f, lb, ub, N, G, s, cv, par, verbose)`: ordenação não dominada, *crowding*, torneio binário por ≺_n, SBX (η_c = 15, p_c = 0.9) e mutação polinomial (η_m = 20, p_m = 1/n), elitismo P_t ∪ Q_t; restrições pela regra de Deb (`cv`, violação ≥ 0).
- Auxiliares: `NonDominatedSort.m` — `[fronts, rank] = NonDominatedSort(F, cv)`; `CrowdingDistance.m` — `d = CrowdingDistance(F, idx)`; `Hypervolume2D.m` — `hv = Hypervolume2D(F, r)`.
- Exemplo: `ex06_6_nsga2.m`, que reproduz os números dos decks 6.6 e 6.2 (partes determinísticas exatas; NSGA-II com verificação estatística).
- Com biblioteca, como no slide «No computador»: `ex06_6_gamultiobj.m` (`gamultiobj`, Global Optimization Toolbox) — **não testado no Octave** (o Octave não tem `gamultiobj`; aí o script só avisa). `NSGA2.m` não depende de toolboxes.
- Apontamentos: `notes/pt/` (deck 6.6).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex06_6_nsga2
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex06_6_nsga2"`; demora cerca de 70 s — são 40 corridas do NSGA-II). A ajuda das funções: `help NSGA2`, `help NonDominatedSort`, `help CrowdingDistance`, `help Hypervolume2D`.

## Método estocástico: verificação estatística

`s` é a semente: a função faz `rng(s)` antes de gerar P₀ (com `s = []` usa o estado atual do gerador). Os números dos slides foram gerados em Python com `numpy.random.default_rng`; o gerador do MATLAB/Octave é outro, por isso **as mesmas sementes dão populações diferentes** e os números do NSGA-II não coincidem um a um. O exemplo verifica:

- n_f = N + G·N = 40 + 50 × 40 = 2040 em todas as corridas (não depende do gerador);
- HV final (exemplo convexo, r = (1.1; 1.65)): a mediana de 20 corridas (`rng(0)` … `rng(19)`) está a menos de 0.005 de 1.637, e 1.637 está entre o mínimo e o máximo;
- geração 0 com 6 não dominados admissíveis: 6 está entre o mínimo e o máximo de 200 populações iniciais;
- tabela «DMS vs. NSGA-II» (coluna do NSGA-II, 20 corridas `rng(100)` … `rng(119)`): em cada orçamento, a mediana está a menos de max(0.005, (Q3 − Q1)/2) da mediana do slide.

A reprodução exata, número a número, está no exemplo Python.

## O que o exemplo imprime

1. As oito soluções dos frames «rank» e «crowding»: rank e d_p de cada uma; quem sai quando só cabem 3 dos 4 de F₁.
2. Hipervolumes do deck 6.2 (três pontos; convergência/diversidade; frentes exatas).
3. NSGA-II no exemplo convexo (N = 40, G = 50): uma corrida com o histórico `[t n_f |F_1| inadm HV]` de 5 em 5 gerações, e as 20 corridas.
4. A coluna do NSGA-II da tabela «DMS vs. NSGA-II» (exemplo não convexo, N = 40, 49 gerações = 2000 avaliações), ao lado da do slide. A coluna do DMS é copiada do slide (o DMS está em `06_5_dms`).

## Saída esperada (resumo)

- F₁ = {1,2,4,8}, F₂ = {3,5,6,7}; d₂ = 1.400, d₄ = 1.233, d₁ = d₈ = Inf; ficam {1,2,8}, sai o 4.
- HV: 0.549 (r = (1.1; 1.3)); 0.654, 0.656, 1.095; frente convexa 1.6483 (exata, 1.1·1.65 − 1/6), não convexa 1.1483 (1.1·1.65 − 2/3).
- **Hipervolume** (deck 6.2): «espalhados, mas longe» dá 0.656; só contam os pontos que dominam r (o deck foi corrigido de 0.652 para 0.656).
- NSGA-II convexo: n_f = 2040; em Octave 8.4, HV mediana 1.6356 [1.6349; 1.6359] (mín. 1.6319, máx. 1.6373); geração 0: entre 4 e 15 não dominados (mediana 8).
- Tabela, em Octave 8.4: 0.996, 1.082, 1.125, 1.132, 1.131 (slide: 0.986, 1.083, 1.125, 1.132, 1.132).
- Última linha: `confere com os slides (verificação estatística nas partes 3 e 4): sim`.

**Contagens.** `info.nfev` conta as avaliações do vetor **f** (cada uma dá f₁, …, f_m): N para P₀ e N por geração, n_f = N + G·N; a população inicial conta. As avaliações da violação `cv` fazem-se nos mesmos pontos e não contam à parte. `info.ngev = info.nhev = 0`. `info.history` tem uma linha por população P_t, t = 0..G: `[t, n_f, |F_1|, inadm, HV]`, com |F_1| = não dominados admissíveis e HV da frente F₁ (com `par.ref`; NaN sem ele). Nota: `gamultiobj` e o `pymoo` contam as gerações de outra forma (no `pymoo`, a população inicial é a 1.ª das 50 gerações: 2000 avaliações).

**Regra de Deb** (`NonDominatedSort` com `cv`): os admissíveis ordenam-se em frentes por dominância; cada inadmissível forma uma frente sozinho, depois, por violação crescente. Sem penalização a calibrar.

Testado em GNU Octave 8.4 (exceto `ex06_6_gamultiobj.m`); `NSGA2` e auxiliares só usam funções comuns a MATLAB e Octave (sem Statistics Toolbox: os quartis são calculados no exemplo por interpolação linear, como `numpy.percentile`).
