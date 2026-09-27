# 5.2 — Simulated annealing (MATLAB)

Estado: **disponível**.

- Função: `SimulatedAnnealing.m` — `[x, fx, info] = SimulatedAnnealing(f, x0, T0, c, sigma, Nmax, lb, ub)`, com o núcleo do slide «No computador»: vizinho `x + sigma*randn(size(x))` projetado em X com `min`/`max`, critério de Metropolis `df <= 0 || rand < exp(-df/T)`, melhor visitado, `T = c*T`. Devolve o **melhor ponto visitado**, `info.nfev` e `info.nacc`.
- Exemplo: `ex05_2_simulated_annealing.m`, que verifica os números do deck 5.2 (ver «Verificação estatística»).
- Apontamentos: `notes/pt/` (deck 5.2).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex05_2_simulated_annealing
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex05_2_simulated_annealing"`; cerca de 25 s). A ajuda: `help SimulatedAnnealing`. Para repetir uma corrida, fixar a semente antes da chamada: `rng(s)`. Com `completo = true` (1.ª linha do exemplo), acrescenta a variação de c no Rastrigin (cerca de 50 s).

## Verificação estatística

Os números dos slides vêm do script Python das figuras (gerador `numpy`, sementes 3, 99 e 11). O MATLAB/Octave tem outro gerador: com `rng(s)` as corridas são reprodutíveis, mas **não são as do slide**. Por isso o exemplo:

- confere **exatamente** o que é determinístico: T_k em k = 30, 100, 200 (1.7; 0.4; 0.05), T_Nmax da tabela, n_f = N_max + 1 (401, 3001) e a pesquisa direta final (37 avaliações a partir do x_SA = 2.4642 do slide);
- confere as **taxas de sucesso estatisticamente**. A taxa do slide também é uma estimativa (em R_s corridas), por isso a diferença entre as duas estimativas tem de ficar abaixo de 3 desvios-padrão binomiais, 3·√(p(1−p)(1/R + 1/R_s)), com p(1−p) ≥ 1/R. Com R = R_s = 100 a tolerância é larga (cerca de 20 pontos percentuais para p = 0.35): verifica a ordem de grandeza, não o valor exato;
- mostra a corrida isolada do «passo a passo» com `rng(3)`, mas o seu resultado (escapar ou não da bacia de 4.72) **não entra na verificação**: é outra corrida (em Octave 8.4 fica na bacia de 4.72; no slide escapa em k = 315).

O exemplo Python (`code/python/05_global_optimization/05_2_simulated_annealing/`) usa o mesmo gerador e as mesmas sementes do script e reproduz os valores **exatos**.

## O que o exemplo imprime

1. A execução passo a passo (x₀ = 6, T₀ = 3, c = 0.98, σ = 0.5, N_max = 400): T, x_k, f(x_k), f_best em k = 0, 30, 100, 200, 400; n_f, aceites; e «SA → método local» (pesquisa direta 1D, passo inicial 0.01, até 10⁻⁶).
2. A tabela das 100 corridas de x₀ = 6 (sucesso |x − x*| < 0.3): σ, T₀, c, N_max, T_Nmax, taxa obtida, taxa do slide, e se estão dentro de 3 d.p.
3. Rastrigin 2D: uma corrida de (4, −4) com `rng(3)` e a chamada do slide «No computador» (30 corridas, `rng(s)`, s = 1..30): mediana, quartis (calculados sem `prctile`, que no MATLAB é da Statistics Toolbox) e sucesso f < 1.

## Saída esperada (resumo, Octave 8.4)

- Passo a passo: T(30) = 1.6698, T(100) = 0.4060, T(200) = 0.0538, n_f = 401; com `rng(3)` o Octave termina em x = 4.7213 (bacia local; o slide, com o gerador do Python, escapa). Pesquisa direta: 37 avaliações, |x − x*| = 5.2·10⁻⁷.
- Tabela (100 corridas): 39, 15, 71, 14, 93, 70 % (slide: 35, 19, 74, 20, 88, 73 %); T_Nmax iguais (9.3e-04, 1.5e-18, 4.0e-01, 9.3e-05, 9.3e-04, 8.5e-18).
- Rastrigin: n_f = 3001; chamada do slide: mediana f = 0.032, sucesso 100 % (slide: 100 %).
- Última linha: `confere com os slides (valores determinísticos exatos; taxas: estatisticamente): sim`.

**Contagens.** `info.nfev` = N_max + 1: conta f(x₀) e uma avaliação por iteração; `info.nacc` = propostas aceites; `info.ngev = info.nhev = 0`; `info.nit` = N_max; `info.T` = T₀c^N_max; `info.xk`, `info.fk` = ponto corrente final. `info.history` tem uma linha por k = 0, …, N_max: `[k T aceite f(x') f(x_k) f_best x'(1:n) x_k(1:n)]` (T é a temperatura usada na iteração k; x' a proposta; x_k o ponto corrente depois da iteração k; na linha 0, x' = x₀).
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave (`rng` existe no Octave desde a versão 7).
