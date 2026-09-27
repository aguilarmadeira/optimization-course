# 3.2.1 — Pesquisa aleatória pura e localizada (MATLAB)

Estado: **disponível**.

- Funções: `RandomSearch.m` — `[x, fx, info] = RandomSearch(f, a, b, N, s, ftarget, verbose)`; `LocalRandomSearch.m` — `[x, fx, info] = LocalRandomSearch(f, x0, r0, m, gamma, tolx, kmax, ftarget, s, verbose)`
- Exemplo: `ex03_2_1_random_search.m`, que verifica estatisticamente os números do deck 3.2.1.
- Apontamentos: `notes/pt/` (deck 3.2.1).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex03_2_1_random_search
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex03_2_1_random_search"`; demora cerca de 12 s). A ajuda das funções: `help RandomSearch`, `help LocalRandomSearch`.

## Método estocástico: verificação estatística

`s` é a semente: a função faz `rng(s)` antes de gerar os pontos (com `s = []` usa o estado atual do gerador). Os números dos slides foram gerados em Python com `numpy.random.default_rng`; o gerador do MATLAB/Octave é outro, por isso **as mesmas sementes dão amostras diferentes** e os números não coincidem um a um. O exemplo faz 30 corridas (`rng(0)` … `rng(29)`), como o slide, e verifica:

- pesquisa pura: a mediana de f está a menos de um fator 3 da mediana do slide (a mediana de 30 corridas varia cerca de 25 %);
- localizada, custo até f < 10⁻⁴: a mediana está a menos de 15 % de 975;
- localizada, melhor f com n_f = 3001 (cerca de 10⁻¹²): a mediana tem a mesma ordem de grandeza (fator 10) da do slide;
- as corridas isoladas do slide (f = 0.053 com N = 1000; f = 5.9·10⁻¹³ na localizada) estão entre o mínimo e o máximo das 30 corridas.

A corrida do slide com N = 100 (f = 1.406) não entra na verificação: com 100 amostras, P(f_best > 1.406) ≈ 1 %. A reprodução exata, número a número, está no exemplo Python.

## O que o exemplo imprime

1. ‖x‖² em [−5,5]², uma corrida (`rng(1)`), N = 1000: as melhorias do melhor ponto (`info.history`).
2. A tabela da análise experimental (N = 10, 100, 1000, 10 000; 30 corridas), ao lado da do slide.
3. A localizada no Rosenbrock de (−1.5, 2), m = 20, r₀ = 1, γ = 0.9: uma corrida de 150 iterações (n_f = 3001) e o custo até f < 10⁻⁴ em 30 corridas.
4. Pura vs. localizada com o mesmo orçamento (3000), 30 corridas.

## Saída esperada (resumo)

- n_f = N; localizada com 150 iterações: n_f = 1 + 20 × 150 = 3001 (inclui f(x₀)); chega a (1, 1).
- Custo até f < 10⁻⁴ (`info.nhit`): mediana perto de 975 (slide: 975 [867; 1052]); em Octave 8.4 dá 1030 [850; 1217].
- Última linha: `confere com os slides (verificação estatística): sim`.

**Contagens.** `info.nfev` conta todas as avaliações de f. Na localizada inclui f(x₀): `nfev = 1 + m·nit`. `info.nhit` é a primeira avaliação com f < `ftarget` (f(x₀) é a 1.ª; as m amostras de cada iteração contam pela ordem em que são avaliadas; `Inf` se nunca). É o *first hitting time* da comparação final do capítulo 3. `info.history`: melhorias do melhor ponto (pura) ou uma linha por iteração `[k, n_f, r, f(x), x]` (localizada).
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave (sem Statistics Toolbox: os quartis são calculados no exemplo por interpolação linear, como `numpy.percentile`).
