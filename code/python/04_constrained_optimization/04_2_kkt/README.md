# 4.2 — Condições KKT (Python)

Estado: **disponível**.

- Não há função da UC neste deck: o código é um exemplo, com uma função auxiliar pequena.
- Auxiliar: `kkt_check.py` — `ok, res = kkt_check(df, g, Jg, u, h=None, Jh=None, lam=None, tol=1e-8)`: verifica as condições KKT num ponto, na convenção da UC (g_j ≥ 0, h_ℓ = 0, 𝓛 = f − uᵀg − λᵀh). `res` tem os resíduos `estac`, `admis`, `compl`, `sinal` e o conjunto ativo `ativas` (índices a começar em 1, como nos slides). Só usa `numpy`; não verifica a LICQ nem classifica o ponto.
- Exemplo: `ex04_2_kkt.py`, que reproduz os números do deck 4.2 (os mesmos que o exemplo MATLAB):
  exemplo 2: u₁ = 1/3, u₂ = 2/3.
- Solvers: `scipy.optimize.minimize(method='trust-constr')`, multiplicadores em `r.v`; `scipy.optimize.linprog` para a produção (como no 1.1).
- Apontamentos: `notes/pt/` (deck 4.2).

## Como correr

Precisa de `numpy` e `scipy`. Na pasta desta secção:

```bash
python ex04_2_kkt.py
```

## O que o exemplo faz

**(a)** Os mesmos oito blocos que o exemplo MATLAB (ver o README MATLAB): círculos largo/apertado e o contra-exemplo u < 0; os 4 casos do exemplo 1; ativa com u = 0; o exemplo 2 (u₁ = 1/3, u₂ = 2/3); os três pontos KKT de min −x²; a sensibilidade (f*(3,1) = 0,934) e os preços-sombra da produção (10, 20, 0); a falha da LICQ; os exercícios 1, 2 e 4.

**(b)** `trust-constr` com `NonlinearConstraint(g, 0, np.inf)` (e `NonlinearConstraint(h, 0, 0)` quando há igualdades), gradientes analíticos (`jac=`), `gtol = 1e-10`, `xtol = 1e-12`. O SciPy escreve 𝓛 = f + vᵀc, logo:

- **u_UC = −`r.v[0]`** (desigualdades g ≥ 0) e **λ_UC = −`r.v[1]`** (igualdades);
- no `fmincon` seria u_UC = `lambda.ineqnonlin` e λ_UC = −`lambda.eqnonlin` (deck 4.2, «No computador»).

A produção com `linprog(c=-c, ...)`: `res.ineqlin.marginals` = (−10, −20, 0) e os preços-sombra são o simétrico (o SciPy minimiza −L), como diz o 1.1.

O `trust-constr` é um método de pontos interiores: no exemplo 2 fica a ~2·10⁻⁵ da solução exata (os multiplicadores das restrições inativas saem ~10⁻⁶, não 0). Por isso o solver é comparado com tolerância 10⁻⁴; os valores com casas decimais nos slides comparam-se à casa mostrada. O aviso «delta_grad == 0.0» (restrições lineares) é silenciado: é esperado.

## Saída esperada (resumo)

- Exemplo 1: tabela dos 4 casos; sobrevive o caso 3, x* = (5/3, 1/3), u₁ = 2/3, f* = 1.
- Exemplo 2: u = (0.3333, 0.6667, 0, 0), f = 2.
- Sensibilidade: f*(3,1) = 0.93444; produção: u = (10, 20, 0), L* = 2600.
- Exercício 4: r* = 4.1841, A* = 267.7384, u* = 5.1885.
- `trust-constr` e `linprog`: os mesmos valores (com u = −r.v, preços-sombra = −marginals).
- Última linha: `confere com os slides: sim`.

Testado com Python 3.11, numpy 2.4 e SciPy 1.17.
