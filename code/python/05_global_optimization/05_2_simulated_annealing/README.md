# 5.2 — Simulated annealing (Python)

Estado: **disponível**.

- Módulo/função: `simulated_annealing.py` — `simulated_annealing(f, x0, T0, c, sigma, Nmax, lb=-inf, ub=inf, rng=None)` → `SAResult`, com o núcleo do slide «No computador» (vizinho gaussiano projetado em X, Metropolis, melhor visitado, fazer T = cT). Devolve o **melhor ponto visitado**, `nfev` e `nacc`.
- `rng`: um `np.random.Generator` (p. ex. `np.random.default_rng(3)`) ou uma semente inteira; o mesmo gerador pode ser passado a várias corridas seguidas, como nos scripts das figuras.
- Exemplo: `ex05_2_simulated_annealing.py`, que reproduz **exatamente** os números do deck 5.2.
- Apontamentos: `notes/pt/` (deck 5.2).

## Como correr

Só precisa de `numpy`. Na pasta desta secção (cerca de 10 s):

```bash
python ex05_2_simulated_annealing.py
```

## Reprodução exata

Os números dos slides vêm de `make_figs_5.py` (gerador `numpy.random.default_rng`, sementes 3, 99 e 11). A função consome os números aleatórios pela mesma ordem do script (uma normal por coordenada para o vizinho; um uniforme só quando Δf > 0, porque o `or` não avalia o segundo termo) e a função é avaliada da mesma maneira, por isso os resultados coincidem bit a bit. O exemplo MATLAB/Octave, com outro gerador, faz uma verificação estatística (ver o README da pasta MATLAB).

## O que o exemplo imprime

1. A execução passo a passo (x₀ = 6, T₀ = 3, c = 0.98, σ = 0.5, N_max = 400, semente 3): T, x_k, f(x_k), f_best em k = 0, 30, 100, 200, 314, 315, 400; o último salto aceite; e «SA → método local» (pesquisa direta 1D, passo inicial 0.01, até 10⁻⁶).
2. A tabela das 100 corridas de x₀ = 6 (semente 99, sucesso |x − x*| < 0.3).
3. Rastrigin 2D: a corrida de (4, −4) (semente 3); as 30 corridas de pontos aleatórios (semente 11) para c = 0.9, 0.95, 0.99, 0.995, 0.999; e, a título informativo, a chamada do slide «No computador» com `default_rng(s)`, s = 1..30.

## Saída esperada (resumo)

- Passo a passo: T(30) = 1.6698, T(100) = 0.4060 (x_k = 4.797, na bacia de 4.72, f_best = 0.893), T(200) = 0.0538; último salto aceite em k = 315, de 4.73 para 3.07 (3.3σ); x = 2.4642 (x* = 2.4681), n_f = 401, 124 aceites. Pesquisa direta: 37 avaliações, |x − x*| = 5.2·10⁻⁷.
- Tabela: 35, 19, 74, 20, 88, 73 % e T_Nmax = 9.3e-04, 1.5e-18, 4.0e-01, 9.3e-05, 9.3e-04, 8.5e-18 — iguais às do slide.
- Rastrigin: de (4, −4), f = 0.0677, |x| = 0.0185, n_f = 3001; sucesso 97, 97, 97, 100, 90 % — iguais aos do slide. Chamada do slide (informativo): mediana f = 0.028, sucesso 97 %.
- Última linha: `confere com os slides: sim`.

O resultado tem os campos `x` (melhor visitado), `fx`, `nit` (= N_max), `nfev` (= N_max + 1, conta f(x₀)), `nacc`, `T` (final), `xk`, `fk` (ponto corrente final), `ngev = nhev = 0`, `history`, `cols`, `flag` e `message`, com as mesmas contagens que a versão MATLAB. `history` tem uma linha por k = 0, …, N_max: `[k, T, aceite, f(x'), f(x_k), f_best, x'_1..x'_n, x_k,1..x_k,n]`.
Testado com Python 3.11 e numpy 2.4.
