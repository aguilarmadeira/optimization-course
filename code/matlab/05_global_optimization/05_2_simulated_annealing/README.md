# 5.2 — Simulated annealing (MATLAB)

Estado: **disponível**.

- Função: `SimulatedAnnealing.m` — `[x, fx, info] = SimulatedAnnealing(f, x0, T0, alpha, n, Ts, sigma, lb, ub, nfmax)`, com o algoritmo das aulas (Deb, 2012, cap. 6). É o núcleo do slide «No computador»:
  - vizinho `x + sigma*randn(size(x))`, projetado em X com `min`/`max`;
  - critério de Metropolis `df < 0 || rand < exp(-df/T)`;
  - `t = t + 1` só quando se aceita;
  - termina quando `T < Ts`; senão, se `mod(t, n) == 0`, faz `T = alpha*T`.

  Devolve o **ponto final** (como nas aulas) e também o **melhor visitado** (`info.xbest`), `info.nfev` e `info.nacc`. `lb`, `ub` e `nfmax` são opcionais (−Inf, Inf e 100 000).
- Exemplo: `ex05_2_simulated_annealing.m`, que verifica os números do deck 5.2 (ver «Verificação estatística»).
- Apontamentos: `notes/pt/` (deck 5.2).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex05_2_simulated_annealing
```

Em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex05_2_simulated_annealing"` (cerca de um minuto).

- A ajuda: `help SimulatedAnnealing`.
- Para repetir uma corrida, fixar a semente antes da chamada: `rng(s)`.
- Com `completo = true` (1.ª linha do exemplo), acrescenta α = 0,3 e 0,8 no Rastrigin (mais cerca de um minuto).

## Verificação estatística

Os números dos slides vêm do script Python das figuras (gerador `numpy`, sementes 3, 99 e 11). O MATLAB/Octave tem outro gerador: com `rng(s)` as corridas são reprodutíveis, mas **não são as do slide**. Por isso o exemplo:

- **Confere exatamente o que é determinístico:**
  - o número de pontos aceites e a temperatura final (191 e 3·0,8¹⁹ = 0,043 na função 1D; 81 e 20·0,5⁸ = 0,078 no Rastrigin);
  - a pesquisa direta final (38 avaliações a partir do x_t = 2,4574 do slide).
- **Confere as taxas de sucesso estatisticamente.**
  - A taxa do slide também é uma estimativa (em R_s corridas). Por isso a diferença entre as duas estimativas tem de ficar abaixo de 3 desvios-padrão binomiais, 3·√(p(1−p)(1/R + 1/R_s)), com p(1−p) ≥ 1/R.
  - Com R = R_s = 100 a tolerância é larga (cerca de 20 pontos percentuais para p = 0,6): verifica a ordem de grandeza, não o valor exato.
  - A mediana de n_f tem de cair no intervalo interquartil [Q1; Q3] do slide.
- **Não verifica a corrida isolada do «passo a passo».** Mostra-a com `rng(3)`, mas o seu resultado (escapar ou não da bacia de 4,72) é outro: em Octave 8.4 fica na bacia de 4,72; no slide escapa em n_f = 325.

O exemplo Python (`code/python/05_global_optimization/05_2_simulated_annealing/`) usa o mesmo gerador e as mesmas sementes do script e reproduz os valores **exatos**.

## O que o exemplo imprime

1. **A execução passo a passo** (x₀ = 6, T₀ = 3, α = 0,8, n = 10, T_s = 0,05, σ = 0,5):
   - o ponto final, o melhor visitado, n_f, os aceites e T final;
   - «SA → método local»: pesquisa direta 1D, com passo inicial 0,01, até 10⁻⁶.
2. **A tabela das 100 corridas** de x₀ = 6:
   - sucesso |x − x*| < 0,3 do ponto final e do melhor visitado;
   - a mediana de n_f;
   - os valores do slide e se estão dentro de 3 d.p.
3. **Rastrigin 2D:**
   - uma corrida de (4, −4) com `rng(3)`;
   - a chamada do slide «No computador» (30 corridas, `rng(s)`, s = 1..30): mediana de f, quartis (calculados sem `prctile`, que no MATLAB é da Statistics Toolbox), sucesso f < 1 e mediana de n_f.

## Saída esperada (resumo, Octave 8.4)

- **Passo a passo:**
  - com `rng(3)`, o Octave termina em x_t = 4,6172 (bacia local; o slide, com o gerador do Python, escapa), com n_f = 501, 191 aceites e T = 0,0432;
  - pesquisa direta: 38 avaliações, |x − x*| = 2,1·10⁻⁷.
- **Tabela (100 corridas):**
  - sucesso (ponto final / melhor visitado) de 41/45, 13/13, 71/95, 16/16, 66/90 e 63/94 %; no slide, 60/63, 11/11, 69/93, 15/15, 60/92 e 64/93 %;
  - n_f mediano de 586, 148, 3318,5, 400, 1484,5 e 3960; no slide, 693, 148, 3318, 404, 1518 e 3942.
- **Rastrigin:** de (4, −4), 81 aceites; na chamada do slide, sucesso de 93 % (slide: 90 %) e n_f mediano de 4912.
- Última linha: `confere com os slides (valores determinísticos exatos; taxas: estatisticamente): sim`.

**Contagens.**

- `info.nfev`: conta f(x₀) e uma avaliação por proposta.
- `info.nacc = info.nit`: pontos aceites (t final).
- `info.ngev = info.nhev = 0`.
- `info.T`: temperatura final; `info.nT`: número de reduções.
- `info.xbest` e `info.fbest`: melhor visitado.
- `info.flag`: 0 se T < T_s; 1 se atingiu `nfmax`.

`info.history` tem uma linha por avaliação: `[n_f t T aceite f(x') f(x_t) f_best x'(1:n) x_t(1:n)]`. A linha 1 é o ponto inicial.

Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave (`rng` existe no Octave desde a versão 7).
