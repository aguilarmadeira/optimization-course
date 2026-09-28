# 6.4 — Método do ε-constrangimento (Python)

Estado: **disponível**.

- Módulo/função: `epsilon_constraint.py` — `epsilon_constraint(f1, f2, x0, bounds, E, verbose=False, grafico=False)` → `EpsResult`. Para cada `e` de `E` resolve min f₂(x) s.a. g(x) = e − f₁(x) ≥ 0 com `scipy.optimize.minimize(method='SLSQP')` (como no slide; a UC não reimplementa o SLSQP) e arranque a quente.
- `EpsResult`: `x`, `fx`, `u` (multiplicadores de f₁ ≤ e: `r.multipliers` do SLSQP nas versões recentes do SciPy; senão, `r.v` do `trust-constr`, como no slide), `nit`, `nfev` (soma de `r.nfev` = avaliações de f₂), `nfev1` (chamadas a f₁), `ncalls` (= nfev + nfev1), `ngev = nhev = 0`, `history` (`[e, x, f1, f2, u, nf2, nf1]`), `cols`, `status`, `flag`, `message`.
- Exemplo: `ex06_4_epsilon_constraint.py`, que reproduz os números do deck 6.4 (os mesmos que o exemplo MATLAB, mais as contagens do SLSQP).
- Apontamentos: `notes/pt/` (deck 6.4).

## Como correr

Precisa de `numpy` e `scipy`. Na pasta desta secção:

```bash
python ex06_4_epsilon_constraint.py
```

## O que o exemplo faz

Os mesmos cinco pontos que a versão MATLAB (à mão; varrimento de 10 valores de e; taxa de troca; restrição inativa; Nantes–Lille). No varrimento, `e = np.linspace(0.05, 0.95, 10)` e x₀ = (0.5, 0.3), a quente: **60 avaliações de f₂ (6 por e), 90 de f₁ (9 por e), 150 chamadas**, como no slide. Mostra também u = `r.v[0]` com `trust-constr` em e = 0.5 (0.9997).

**Forma da restrição e contagens.** A restrição escreve-se na forma da UC, `{'type': 'ineq', 'fun': lambda x: e - f1(x)}` (g ≥ 0), como no slide e no script das figuras: dá 60 avaliações de f₂ e 90 de f₁ (150 chamadas).

## Saída esperada (resumo)

- Tabela do varrimento: x₁ = e, x₂ = 0, f₂ = 1 − e², u = 2e, `nf2 = 6`, `nf1 = 9` em cada linha.
- `60 avaliações de f2 (soma de r.nfev), 90 de f1 (150 chamadas)`.
- `u = 1.0000`, `f2* = 0.7399; diferença -0.0101`, `u = 1.8000 … compra 0.0181`.
- Última linha: `confere com os slides: sim`.

As contagens 60/90/150 são as do SLSQP do SciPy e só se verificam em Python; o MATLAB/Octave usa outro solver e dá os mesmos pontos e multiplicadores com outras contagens.
Testado com Python 3.11, numpy 2.4 e SciPy 1.17.
