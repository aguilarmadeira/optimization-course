# 4.2 — Condições KKT (MATLAB)

Estado: **disponível**.

- Não há função da UC neste deck: o código é um exemplo, com uma função auxiliar pequena.
- Auxiliar: `KKTCheck.m` — `[ok, res] = KKTCheck(df, g, Jg, u, h, Jh, lam, tol)`: verifica as condições KKT num ponto, na convenção da UC (g_j ≥ 0, h_ℓ = 0, 𝓛 = f − uᵀg − λᵀh). Recebe o gradiente de f, os valores e as jacobianas das restrições (linha j = ∇g_jᵀ) e os multiplicadores; devolve em `res` os resíduos `estac` (‖∇f − J_gᵀu − J_hᵀλ‖∞), `admis` (max(−g_j, |h_ℓ|)), `compl` (max |u_j g_j|), `sinal` (max(−u_j, 0)) e o conjunto ativo `ativas`. Não verifica a LICQ nem classifica o ponto. `help KKTCheck`.
- Exemplo: `ex04_2_kkt.m`, que reproduz os números do deck 4.2:
  exemplo 2: u₁ = 1/3, u₂ = 2/3.
- Solvers: `fmincon` (`lambda.ineqnonlin`, `lambda.eqnonlin`) e `linprog` (`lambda.ineqlin`), Optimization Toolbox, só se existirem; em GNU Octave, `sqp` e `glpk` (Octave base).
- Apontamentos: `notes/pt/` (deck 4.2).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex04_2_kkt
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex04_2_kkt"`).

## O que o exemplo faz

**(a) Pelas condições KKT**, com as contas do slide, e `KKTCheck` em cada ponto:

1. «De 4.1 para 4.2»: círculo largo (x* = (1,1), inativa, u* = 0) e apertado (u* = √2 − 1); «Geometria (2)»: (√2, √2) no círculo largo tem u = 1/√2 − 1 < 0 e **não** é ponto KKT (falha o sinal).
2. Exemplo 1: os 2^J = 4 casos da receita, cada um um sistema linear em (x₁, x₂, u₁, u₂); a tabela do slide e o caso 3 sobrevivente, x* = (5/3, 1/3), u₁ = 2/3, f* = 1; ∇f(x*) = (2/3)(−1,−4).
3. «Ativa não implica u > 0»: min x² s.a. x ≥ 0.
4. Exemplo 2: em (2,1) as ativas são {1, 2}; resolve ∇f = u₁∇g₁ + u₂∇g₂ → u = (1/3, 2/3, 0, 0), LICQ (característica 2), f = 2.
5. «KKT dá candidatos»: min −x² s.a. x + 1 ≥ 0, 2 − x ≥ 0 — três pontos KKT (x = 0 máximo, −1 mínimo local, 2 mínimo global).
6. «O que medem os multiplicadores»: exemplo 1 com b = 3,1 (f* = 0,934; aproximação 1 − (2/3)·0,1 = 0,933) e a produção do 1.1 em g ≥ 0 (preços-sombra u = (10, 20, 0) pelas duas ativas: 2u₁ + u₂ = 40, u₁ + u₂ = 30).
7. «Quando KKT falha: LICQ»: em (1,0) os gradientes ativos são dependentes; resíduo mínimo 1, não há u.
8. Exercícios 1, 2 e 4 (as respostas a cinzento): ex. 1 x* = (0,3), ativas {1,2} (u = (1,1,0), não está no slide); ex. 2 x* = (4/5, 8/5), λ = 8/5, desigualdades inativas; ex. 4 (lata com h ≤ 6) r* = 4,18, A* = 267,7, u* ≈ 5,19 (λ = 0,5723, não está no slide).

**(b) Com um solver** (exemplo 1, exemplo 2, exercícios 2 e 4; produção em PL):

| Solver | Chamada | Conversão para a UC |
|---|---|---|
| `fmincon` (MATLAB) | `nonlcon = @(x) deal(-g(x), h(x))` (c = −g ≤ 0, c_eq = h) | **u = `lambda.ineqnonlin`**; **λ = −`lambda.eqnonlin`** (𝓛_fmincon = f + μᵀc + λ_eqᵀc_eq) |
| `sqp` (Octave) | `sqp(x0, f, h, g)` com g(x) ≥ 0, h(x) = 0 | nenhuma: `lambda = [λ; u]` já na convenção da UC (verificado nos exemplos) |
| `linprog` (MATLAB) | `linprog(-c, A, b, [], [], lb)` (min −L) | preços-sombra = `lambda.ineqlin` = (10, 20, 0) |
| `glpk` (Octave) | `glpk(c, A, b, lb, [], 'UUU', 'CC', -1)` (maximiza L) | preços-sombra = `extra.lambda` = (10, 20, 0) |

Os códigos de saída do `sqp` 101 e 104 são paragem normal nestes problemas.

## Saída esperada (resumo)

- Exemplo 1: casos 1 (x = (2,1), g₁ = −3), 2 (x = (4/3, 4/3), u₂ = −4/3), 3 (x = (5/3, 1/3), u₁ = 2/3 ✓), 4 (x = (3/5, 3/5), u = (22/25, −48/25)).
- Exemplo 2: u = (0.3333, 0.6667, 0, 0), f = 2.
- Sensibilidade: f*(3,1) = 0.93444; produção: u = (10, 20, 0), L* = 2600.
- Exercício 4: r* = 4.1841, A* = 267.7384, u* = 5.1885.
- Solver (Octave `sqp`/`glpk`): os mesmos valores.
- Última linha: `confere com os slides: sim`.

Nota: no caso 2 do exemplo 1 falham duas condições (g₁ = −11/3 < 0 e u₂ = −4/3 < 0); o slide só aponta u₂ < 0 — basta uma para o excluir.

Testado em GNU Octave 8.4 (ramos `sqp` e `glpk`); os ramos `fmincon`/`linprog` usam só a sintaxe documentada mas não foram corridos aqui (sem MATLAB). Só usa funções comuns a MATLAB e Octave.
