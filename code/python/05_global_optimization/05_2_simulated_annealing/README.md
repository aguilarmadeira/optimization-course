# 5.2 — Simulated annealing (Python)

Estado: **disponível**.

- Módulo/função: `simulated_annealing.py` — `simulated_annealing(f, x0, T0, alpha, n, Ts, sigma, lb=-inf, ub=inf, rng=None, nfmax=100000)` → `SAResult`, com o algoritmo das aulas (Deb, 2012, cap. 6), o do slide «Algoritmo».
- Devolve o **ponto final** `x` (como nas aulas) e também o **melhor visitado** `xbest`, que na prática convém guardar. Devolve ainda `nfev` e `nacc`.
- `rng`: um `np.random.Generator` (p. ex. `np.random.default_rng(3)`) ou uma semente inteira; o mesmo gerador pode ser passado a várias corridas seguidas, como nos scripts das figuras.
- Exemplo: `ex05_2_simulated_annealing.py`, que reproduz **exatamente** os números do deck 5.2.
- Apontamentos: `notes/pt/` (deck 5.2).

## O algoritmo (como nas aulas)

1. **Passo 1.** Escolher x₀, a temperatura inicial T₀ (alta), o número n de pontos aceites a cada temperatura, α (nas aulas entre 0,5 e 0,99) e a temperatura mínima T_s; fazer t = 0.
2. **Passo 2.** Gerar um vizinho aleatório x′ de x_t; aqui, x′ = x_t + σζ, com ζ ~ N(0, I), projetado em X.
3. **Passo 3.** Δf = f(x′) − f(x_t).
   - Se Δf < 0, aceitar: x_{t+1} = x′ e t = t + 1.
   - Senão, gerar r ~ U(0, 1). Se r < e^{−Δf/T}, aceitar; senão, voltar ao passo 2.
4. **Passo 4.** Se T < T_s, terminar. Senão, se t mod n = 0, fazer T = αT. Voltar ao passo 2.

O contador t só avança quando se aceita. A T baixa quase todas as propostas são rejeitadas, e cada patamar custa cada vez mais avaliações. Por isso há também um limite `nfmax` de avaliações (`flag = 1`).

## Como correr

Só precisa de `numpy`. Na pasta desta secção (cerca de 30 s):

```bash
python ex05_2_simulated_annealing.py
```

## Reprodução exata

Os números dos slides vêm de `make_figs_5.py` (gerador `numpy.random.default_rng`, sementes 3, 99 e 11). A função consome os números aleatórios pela mesma ordem do script:

- uma normal por coordenada para o vizinho;
- um uniforme só quando Δf ≥ 0, porque o `or` não avalia o segundo termo.

A função é avaliada da mesma maneira, por isso os resultados coincidem bit a bit. O exemplo MATLAB/Octave usa outro gerador e faz uma verificação estatística (ver o README da pasta MATLAB).

## O que o exemplo imprime

1. **A execução passo a passo** (x₀ = 6, T₀ = 3, α = 0,8, n = 10, T_s = 0,05, σ = 0,5, semente 3):
   - t, T, x_t, f(x_t) e f_best em n_f = 1, 60, 150, 300, 324, 325 e no fim;
   - o salto aceite através da barreira;
   - «SA → método local»: pesquisa direta 1D a partir do ponto final, com passo inicial 0,01, até 10⁻⁶.
2. **A tabela das 100 corridas** de x₀ = 6 (semente 99, T_s = 0,05):
   - sucesso |x − x*| < 0,3 do ponto final e do melhor visitado;
   - mediana e quartis de n_f.
3. **Rastrigin 2D** (T₀ = 20, α = 0,5, n = 10, T_s = 0,1, σ = 0,5):
   - a corrida de (4, −4), com semente 3;
   - as 30 corridas de pontos aleatórios (semente 11) para α = 0,3, 0,5 e 0,8;
   - a título informativo, a chamada do slide «No computador» com `default_rng(s)`, s = 1..30.

## Saída esperada (resumo)

- **Passo a passo:**
  - em n_f = 60, 150 e 300, T = 1,229, 0,322 e 0,106, com o ponto na bacia de 4,72;
  - salto em n_f = 325, de 4,59 para 2,92 (3,3σ), com T = 0,084;
  - termina com T = 0,043 < T_s em x_t = 2,4574 (melhor visitado 2,4692), com n_f = 507 e 191 aceites;
  - pesquisa direta: 38 avaliações, |x − x*| = 2,1·10⁻⁷.
- **Tabela:** sucesso do ponto final / melhor visitado de 60/63, 11/11, 69/93, 15/15, 60/92 e 64/93 %, com n_f mediano de 693, 148, 3318, 404, 1518 e 3942. São os valores do slide.
- **Rastrigin:**
  - de (4, −4): f = 0,183, n_f = 5421, 81 aceites;
  - α = 0,3, 0,5 e 0,8: 97/97, 90/100 e 97/100 %, com n_f mediano de cerca de 3630, 4814 e 18 572. São os valores do slide.
- Última linha: `confere com os slides: sim`.

O resultado tem os mesmos campos e contagens que a versão MATLAB:

- `x` e `fx`: ponto final;
- `xbest` e `fbest`: melhor visitado;
- `nit = nacc = t`;
- `nfev`: conta f(x₀);
- `T` (final) e `nT` (número de reduções);
- `ngev = nhev = 0`;
- `history`, `cols`, `flag` (0: T < T_s; 1: `nfmax`) e `message`.

`history` tem uma linha por avaliação: `[n_f, t, T, aceite, f(x'), f(x_t), f_best, x'_1..x'_n, x_t,1..x_t,n]`.

Testado com Python 3.11 e numpy 2.4.
