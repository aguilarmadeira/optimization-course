# 6.6 — NSGA-II (Python)

Estado: **disponível**.

- Módulo/função: `nsga2.py` — `nsga2(f, lb, ub, N, G, s=None, cv=None, pc=0.9, eta_c=15, pm=None, eta_m=20, ref=None, verbose=False)` → `NSGA2Result`: ordenação não dominada, *crowding*, torneio binário por ≺_n, SBX e mutação polinomial, elitismo P_t ∪ Q_t; restrições pela regra de Deb (`cv`, violação ≥ 0).
- Auxiliares: `non_dominated_sort.py` — `non_dominated_sort(F, cv=None)` → `(fronts, rank)`; `crowding_distance.py` — `crowding_distance(F, idx)`; `hypervolume_2d.py` — `hypervolume_2d(F, r)`.
- Exemplo: `ex06_6_nsga2.py`, que reproduz **exatamente** os números dos decks 6.6 e 6.2 (mesmas sementes e mesmo gerador do script das figuras).
- Alternativa com biblioteca, como no slide «No computador»: `nsga2_exemplo_uc.py` (`pymoo`, `pip install pymoo matplotlib`): 40 pontos, 2000 avaliações, HV 1,636 (no `pymoo` a população inicial conta como 1.ª das 50 gerações; aqui: inicial + 50 gerações = 2040 avaliações, HV 1.637).
- Apontamentos: `notes/pt/` (deck 6.6).

## Como correr

`ex06_6_nsga2.py` só precisa de `numpy`. Na pasta desta secção:

```bash
python ex06_6_nsga2.py        # cerca de 7 s
python nsga2_exemplo_uc.py    # pymoo + matplotlib (grava nsga2_exemplo_uc.png)
```

## Método estocástico: reprodução exata

`s` é a semente (inteiro) ou um `numpy.random.Generator`. `nsga2` faz as mesmas chamadas ao gerador, pela mesma ordem, que a implementação didática de `make_figs_6.py` (P₀ uniforme; em cada par de filhos: dois torneios, `random() < pc`, SBX, duas mutações). Com as sementes do script — 4 (exemplo convexo), 5 (uma corrida no não convexo), 100–119 (20 corridas da tabela) — reproduz os números dos slides até ao último dígito (HV convexo = 1.6368824197851621).

## O que o exemplo imprime

1. As oito soluções dos frames «rank» e «crowding»: rank e d_p de cada uma; quem sai quando só cabem 3 dos 4 de F₁.
2. Hipervolumes do deck 6.2 (três pontos; convergência/diversidade; frentes exatas).
3. NSGA-II no exemplo convexo (N = 40, G = 50, semente 4), com o histórico `[t n_f |F_1| inadm HV]` de 5 em 5 gerações.
4. A coluna do NSGA-II da tabela «DMS vs. NSGA-II» (exemplo não convexo, N = 40, 49 gerações = 2000 avaliações, 20 corridas), ao lado da do slide. A coluna do DMS é copiada do slide (o DMS está em `06_5_dms`).

## Saída esperada (resumo)

- F₁ = {1,2,4,8}, F₂ = {3,5,6,7}; d₂ = 1.400, d₄ = 1.233, d₁ = d₈ = inf; ficam {1,2,8}, sai o 4.
- HV: 0.549 (r = (1.1; 1.3)); 0.654, 0.656, 1.095; frente convexa 1.6483 (exata, 1.1·1.65 − 1/6), não convexa 1.1483 (1.1·1.65 − 2/3).
- **Hipervolume** (deck 6.2): «espalhados, mas longe» dá 0.656; só contam os pontos que dominam r (o deck foi corrigido de 0.652 para 0.656).
- NSGA-II convexo: geração 0 com 6 não dominados (14 inadmissíveis); geração 50: 40 pontos, n_f = 2040, HV = 1.637.
- Tabela: 0.986 [0.968; 1.006], 1.083 [1.066; 1.096], 1.125 [1.117; 1.129], 1.132 [1.131; 1.133], 1.132 [1.132; 1.133] — iguais às do slide.
- Última linha: `confere com os slides: sim`.

**Contagens.** `nfev` conta as avaliações do vetor **f**: N para P₀ e N por geração, n_f = N + G·N; a população inicial conta. As avaliações de `cv` fazem-se nos mesmos pontos e não contam à parte. `ngev = nhev = 0`. O resultado tem os campos `X`, `F` (frente F₁ de P_G), `nit`, `nfev`, `ngev`, `nhev`, `P`, `FP`, `cvP`, `rank`, `crowd`, `history` (uma linha por P_t, t = 0..G: `[t, n_f, |F_1|, inadm, HV]`), `flag` e `message`, como a versão MATLAB (`[X, F, info]`). Os índices das frentes são base 0; o rank começa em 1, como no slide.
Testado com Python 3.11 e numpy 2.4 (e pymoo 0.6.2 para `nsga2_exemplo_uc.py`).
