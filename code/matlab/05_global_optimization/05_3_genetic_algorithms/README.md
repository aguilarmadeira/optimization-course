# 5.3 — Algoritmos genéticos (MATLAB)

Estado: **disponível**.

- Função: `GeneticAlgorithm.m` — `[xb, fb, info] = GeneticAlgorithm(f, n, lb, ub, N, G, pc, pm, sigma, e, s)`, com o núcleo do slide «No computador» (codificação real): elites, torneio, `cruzarBLX`, mutação gaussiana `Q + M.*sigma.*randn(N,n)`, projeção em X, elitismo, avaliação. O 11.º argumento, `s` (tamanho do torneio), é opcional: por omissão 2, e então são as linhas `i1`, `i2` do slide.
- Auxiliar: `cruzarBLX.m` — `Q = cruzarBLX(pais, pc, beta)`: BLX-β aos pares (1-2, 3-4, …), cada par com probabilidade pc; a ~ U(−β, 1+β) por coordenada; β = 0.5 por omissão.
- Exemplo: `ex05_3_genetic_algorithms.m`, que verifica os números do deck 5.3 (ver «Verificação estatística»).
- Apontamentos: `notes/pt/` (deck 5.3).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex05_3_genetic_algorithms
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex05_3_genetic_algorithms"`; cerca de 25 s). A ajuda: `help GeneticAlgorithm`, `help cruzarBLX`. Para repetir uma corrida: `rng(s)` antes da chamada. Com `completo = true` (1.ª linha do exemplo), a tabela dos parâmetros usa 300 corridas por linha, como o slide (cerca de 1 min).

## Verificação estatística

Os números dos slides vêm dos scripts Python das figuras (gerador `numpy`, sementes 5, 2, 21 e 11). O MATLAB/Octave tem outro gerador: com `rng(s)` as corridas são reprodutíveis, mas **não são as do slide**. Por isso o exemplo:

- confere **exatamente** a geração à mão (parte 1: parte dos números sorteados que o slide mostra e refaz as contas) e as contagens n_f = N(G + 1) (2440; 420, 126, 1260);
- confere as **taxas de sucesso estatisticamente**: a diferença entre a taxa obtida (R corridas) e a do slide (R_s = 300 ou 30 corridas, também uma estimativa) tem de ficar abaixo de 3·√(p(1−p)(1/R + 1/R_s)), com p(1−p) ≥ 1/R. Com R = 100 a tolerância é cerca de 15 pontos percentuais (p = 0.75); com R = 300, cerca de 10.6;
- as corridas isoladas (Rastrigin com `rng(2)`) mostram-se, mas o valor de f não entra na verificação.

O exemplo Python (`code/python/05_global_optimization/05_3_genetic_algorithms/`) usa o mesmo gerador e as mesmas sementes dos scripts e reproduz os valores **exatos**.

**Nota sobre a tabela dos parâmetros.** Em 3000 corridas (Python; nas três linhas também testadas em Octave, a diferença foi menor que 1,5 pontos percentuais) as taxas «de longo prazo» são cerca de 75, 61, **54**, **39**, 98 % (referência, s = 6, p_c = 0, N = 6, N = 60), 74 e 78 % (p_m = 0 e 0.8), 97 e **40** % (mesmo orçamento). Os 60 %, 34 % e 46 % do slide são estimativas de 300 corridas a cerca de 2 desvios-padrão desses valores; as conclusões do slide (pesa sobretudo N; com n_f = 420, N = 60 é muito melhor do que N = 6) mantêm-se.

## O que o exemplo imprime

1. A geração à mão (min x² em [−5, 5], 8 bits, N = 6, roleta, p_c = 0.8, p_m = 0.05): tabela da população inicial (cromossoma, x, f, F, p_i, acumulada), pais da roleta, cruzamentos, mutação, x e f dos descendentes, médias 11.9 → 10.4 e o melhor 0.111 → 0.323.
2. Rastrigin 2D (N = 40, G = 60, p_c = 0.9, p_m = 0.1, σ = 0.5, 2 elites): melhor e média nas gerações 0, 5, 20, 60 de uma corrida com `rng(2)`; a chamada do slide «No computador» (30 corridas, `rng(s)`): mediana, quartis (sem `prctile`, que no MATLAB é da Statistics Toolbox) e sucesso f < 1.
3. A tabela «O que fazem os parâmetros?» na função de referência (sucesso |x − x*| < 0.3), com as linhas p_m = 0, p_m = 0.8 e as duas de mesmo orçamento (n_f = 420).

## Saída esperada (resumo, Octave 8.4)

- Geração à mão: todos os valores iguais aos do slide.
- Rastrigin: n_f = 2440; com `rng(2)`, f = 5·10⁻¹⁴; chamada do slide: mediana f = 7.7·10⁻¹², sucesso 100 % (slide: 100 %).
- Tabela (100 corridas, `completo = false`): 67, 62, 54, 37, 98 %; 73, 79 %; 95, 46 % — todas dentro da tolerância; n_f = 420, 126, 1260 iguais. Com `completo = true` (300 corridas): 73.0, 63.3, 50.3, 43.7, 97.0; 73.0, 80.0; 96.7, 43.0 %.
- Última linha: `confere com os slides (valores determinísticos exatos; taxas: estatisticamente, 100 corridas): sim`.

**Contagens.** `info.nfev` = N(G + 1): a população inicial e N avaliações por geração (as elites são reavaliadas, como no slide). `info.ngev = info.nhev = 0`; `info.nit` = G; `info.P`, `info.F` = população final e os seus valores de f. `info.history` tem uma linha por geração t = 0, …, G: `[t  melhor f  média de f  x_melhor(1:n)]`. Como no slide, `F` guarda os valores de f (menor = melhor), não uma aptidão a maximizar. O resultado é o melhor da população final (com elitismo, é também o melhor de todas as gerações).
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave (`rng` existe no Octave desde a versão 7).
