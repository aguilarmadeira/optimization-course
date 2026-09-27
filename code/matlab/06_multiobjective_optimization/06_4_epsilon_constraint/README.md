# 6.4 — Método do ε-constrangimento (MATLAB)

Estado: **disponível**.

- Função: `EpsilonConstraint.m` — `[x, fx, info] = EpsilonConstraint(f1, f2, x0, lb, ub, E, verbose, grafico)`. Para cada `e` de `E` resolve min f₂(x) s.a. g(x) = e − f₁(x) ≥ 0, lb ≤ x ≤ ub, com arranque a quente. O núcleo é o do slide «Em MATLAB: varrer ε com fmincon»: `nonl = @(x) deal(f1(x) - e, [])`, `u = lambda.ineqnonlin`, `info.nfev = info.nfev + out.funcCount`, `x0 = x`.
- `info`: `u` (multiplicador de f₁ ≤ e, a taxa de troca), `nfev` (soma de `out.funcCount` com `fmincon`; chamadas a f₂ com `sqp`), `nf2`, `nf1` e `ncalls` (chamadas a f₂ e a f₁ contadas diretamente, e a soma), `ngev = nhev = 0`, `nit`, `history` (uma linha por e: `[e x f1 f2 u nf2 nf1]`), `cols`, `exitflag`, `solver`, `flag`, `message`. `help EpsilonConstraint`.
- Exemplo: `ex06_4_epsilon_constraint.m`, que reproduz os números do deck 6.4.
- Solver: `fmincon` (Optimization Toolbox), como no slide, se existir; em GNU Octave, `sqp` (Octave base) com a restrição já na forma da UC; o multiplicador é `lambda(1)` (o `sqp` devolve primeiro o da restrição, depois os dos limites). Os pontos e os multiplicadores coincidem.
- Apontamentos: `notes/pt/` (deck 6.4).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex06_4_epsilon_constraint
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex06_4_epsilon_constraint"`).

## O que o exemplo faz

Exemplo não convexo da UC: min (x₁, 1 − x₁² + x₂), x ∈ [0,1] × [0; 0.6], x₀ = (0.5, 0.3).

1. «À mão»: e = 0.25, 0.50, 0.75 → (0.25; 0.9375), (0.50; 0.75), (0.75; 0.4375).
2. Varrimento do slide, `e = 0.05:0.1:0.95` (10 valores; a variável chama-se `e`, não `eps`): 10 pontos x* = (e, 0), f₂ = 1 − e², u = 2e.
3. Taxa de troca: e = 0.5 → u = 1; f₂*(0.51) = 0.7399, diferença −0.0101; e = 0.9 → u = 1.8, a folga até 0.91 compra 0.0181.
4. Restrição inativa: e = 1.05 e 1.2 → a mesma solução (1, 0), u = 0.
5. Nantes–Lille (6.1): o mais barato com duração ≤ e: 4 h → avião (134 €), 5 h → TGV low-cost (40 €), 10 h → autocarro 2 (28 €).

## Contagens

O slide dá, em Python (SLSQP), 10 pontos com **60 avaliações de f₂ e 90 de f₁ (150 chamadas)**. Essas contagens são do SLSQP do SciPy e **só se verificam na versão Python**. Aqui o solver é outro: com o `sqp` do Octave, n_f = 70 chamadas a f₂ e 140 a f₁ (210); com o `fmincon`, a soma de `out.funcCount`. O exemplo imprime-as só a título informativo; não entram na verificação. O ponto inicial não é avaliado à parte; as avaliações finais de f₁ e f₂ para a tabela não contam.

## Saída esperada (resumo)

- Tabela do varrimento: e = 0.05 … 0.95, x₁ = e, x₂ = 0, f₂ = 0.9975, 0.9775, …, 0.0975, u = 0.1, 0.3, …, 1.9.
- `u = 1.0000`, `f2* = 0.7399; diferença -0.0101`, `u = 1.8000 … compra 0.0181`.
- Inativa: `x* = (1.0000; 0.0000) … u = 0.0000`.
- Última linha: `confere com os slides: sim`.

Testado em GNU Octave 8.4 (ramo `sqp`); o ramo `fmincon` usa só a sintaxe do slide mas não foi corrido aqui (sem MATLAB). Só usa funções comuns a MATLAB e Octave.
