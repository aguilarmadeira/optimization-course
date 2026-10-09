# 5.3 — Algoritmos genéticos (Python)

Estado: **disponível**.

- Módulo/funções: `genetic_algorithm.py` — `genetic_algorithm(f, n, lb, ub, N, G, pc, pm, sigma, e, s=2, rng=None)` → `GAResult`, com o núcleo do slide «No computador» (elites, torneio de tamanho s, BLX-0.5, mutação gaussiana, projeção em X, elitismo); e a auxiliar `cruzar_blx(pais, pc, rng, beta=0.5)` (a `cruzarBLX` do slide).
- `rng`: um `np.random.Generator` (p. ex. `np.random.default_rng(2)`) ou uma semente inteira; o mesmo gerador pode ser passado a várias corridas seguidas, como nos scripts das figuras.
- Exemplo: `ex05_3_genetic_algorithms.py`, que reproduz **exatamente** os números do deck 5.3. Usa `simulated_annealing` da pasta `05_2_simulated_annealing/` (no `sys.path` através de `uc_setup`) (a amostra de 30 pontos do script intercala SA e AG no mesmo gerador).
- Apontamentos: `notes/pt/` (deck 5.3).

## Como correr

Só precisa de `numpy`. Na pasta desta secção (cerca de 20 s):

```bash
python ex05_3_genetic_algorithms.py
```

Com `COMPLETO = False` (no início do exemplo), a parte 3 faz só 30 corridas por linha e mostra as taxas sem as conferir: os valores exatos do slide só saem com `COMPLETO = True` (300 corridas, por omissão).

## Reprodução exata

Os números dos slides vêm de `make_figs_5.py` (geração à mão, semente 5; Rastrigin, sementes 2 e 21) e de `make_figs_54.py` (tabela dos parâmetros, semente 11). `genetic_algorithm` consome os números aleatórios pela mesma ordem de `make_figs_5.py`, e os resultados coincidem bit a bit. Há uma diferença de ordem em `make_figs_54.py`: sorteia os candidatos do torneio como uma matriz N×s (linha = pai), enquanto a função (como o slide, com `i1`, `i2`) os sorteia como s×N (linha = candidato). A distribuição é a mesma; só muda a ordem em que os números saem do gerador. Para reproduzir a tabela bit a bit, o exemplo passa um pequeno adaptador do gerador (`OrdemMakeFigs54`) que faz o sorteio na ordem do script e o devolve transposto. Na função de referência 1D, f é avaliada num vetor de 1 componente (`f1(x)[0]`), como o script, que avalia a população de uma vez: assim os valores coincidem até ao último bit.

A geração à mão (parte 1) é o AG binário do slide (roleta, cruzamento a 1 ponto, mutação por bit), não `genetic_algorithm`; está escrita no exemplo, passo a passo, com a semente 5 do script.

## O que o exemplo imprime

1. A geração à mão (min x² em [−5, 5], 8 bits, N = 6): população inicial (cromossoma, x, f, F, p_i, acumulada), ξ e pais da roleta, cortes, descendentes, a mutação (descendente 2, bit 6), x e f dos descendentes, médias e melhor.
2. Rastrigin 2D (N = 40, G = 60, p_c = 0.9, p_m = 0.1, σ = 0.5, 2 elites): melhor e média nas gerações 0, 5, 20, 60 (semente 2); as 30 corridas da amostra do script (semente 21), SA e AG; e, a título informativo, a chamada do slide «No computador» com `default_rng(s)`, s = 1..30.
3. A tabela «O que fazem os parâmetros?» (300 corridas, semente 11, sucesso |x − x*| < 0.3), com p_m = 0, p_m = 0.8 e as duas linhas de mesmo orçamento (n_f = 420).

## Saída esperada (resumo)

- Geração à mão: cromossomas 11010110, 10001000, …; ξ = 0.679, 0.870, 0.227, 0.895, 0.872, 0.019 → pais 2 4 2 5 4 1; cortes após os bits 3, 4, 7; uma mutação; média de f 11.9 → 10.4, melhor 0.111 → 0.323 — tudo igual ao slide.
- Rastrigin: gerações 0/5/20/60: melhor 6.48 / 0.581 / 1.8·10⁻⁵ / 8.2·10⁻¹³, média 31.3 / 15.2 / 2.0 / 1.4 (os valores da figura); f < 10⁻¹², n_f = 2440. Amostra de 30 pontos: AG 100 % (n_f = 2440), SA das aulas 97 % (ponto final; n_f mediana 5323).
- Tabela: 74.7, 60.7, 60.0, 33.7, 96.3 %; 76.3, 78.0 %; 96.3, 45.7 % → 75/61/60/34/96, 76/78, 96/46 % do slide; n_f = 420, 126, 1260.
- Última linha: `confere com os slides: sim`.

O resultado tem os campos `x`, `fx` (o melhor da população final), `nit` (= G), `nfev` (= N(G + 1): população inicial + N por geração, elites reavaliadas), `ngev = nhev = 0`, `P`, `F` (população final e os seus valores de f), `history` (uma linha por geração: `[t, melhor f, média de f, x_melhor_1..n]`), `cols`, `flag` e `message`, com as mesmas contagens que a versão MATLAB.
Testado com Python 3.11 e numpy 2.4.
