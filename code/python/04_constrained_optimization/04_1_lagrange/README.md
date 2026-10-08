# 4.1 — Multiplicadores de Lagrange (Python)

Estado: **disponível**.

- Não há função da UC neste deck: o código é um exemplo.
- Exemplo: `ex04_1_lagrange.py`, que reproduz os números do deck 4.1 (os mesmos que o exemplo MATLAB):
  lata: r* = 3,745, h* = 7,49, A* ≈ 264,4 cm², λ* = −0,534 (dA*/dV = −λ* = 0,534).
- Solver: `scipy.optimize.minimize(method='trust-constr')`, multiplicadores em `r.v` (deck 4.2, «No computador»).
- Apontamentos: `notes/pt/` (deck 4.1).

## Como correr

Precisa de `numpy` e `scipy`. Na pasta desta secção:

```bash
python ex04_1_lagrange.py
```

## O que o exemplo faz

**(a) Pelas condições de Lagrange** (convenção das aulas 𝓛 = f + Σ λ_ℓ h_ℓ), com as contas do slide e a verificação numérica de ∇ₓ𝓛 = 0 e h = 0: o exemplo resolvido (com «Geometria (1)», a sensibilidade df*/dδ = −λ* e a Hessiana orlada), os candidatos (±1, ±1) e a sua classificação (2.ª ordem e Hessiana orlada), a lata, as várias restrições e o exemplo em que Lagrange falha (resíduo mínimo 1). Ver o README MATLAB para o detalhe.

**(b) Com `trust-constr`** (exemplo resolvido, lata, várias restrições), com `NonlinearConstraint(h, 0, 0)` e os gradientes analíticos (`jac=`). O SciPy escreve 𝓛 = f + vᵀc, como nas aulas, logo **λ = `r.v[0]`**, sem troca de sinal — como no `fmincon` (λ = `lambda.eqnonlin`). Tolerâncias `gtol = 1e-10`, `xtol = 1e-12`. O aviso «delta_grad == 0.0» (restrições lineares) é silenciado no script: é esperado e não afeta o resultado.

## Saída esperada (resumo)

- Exemplo resolvido: x* = (0.3333, 0.6667), λ* = −1.3333, f* = 0.6667; Hessiana orlada: det = −6.
- Candidatos: (−1,−1) mínimo (λ = 0.5, det orlada = −8), (1,1) máximo (λ = −0.5, det orlada = 8).
- Lata: r* = 3.7449 cm, h* = 7.4899 cm, A* = 264.3568 cm², λ* = −0.5341 cm²/cm³.
- Várias restrições: x* = (1.5, 0.5, 1), λ* = (−2, −1), f* = 3.5.
- `trust-constr`: os mesmos valores, com λ = r.v (a lata precisa de ~136 iterações, paragem por `xtol`).
- Última linha: `confere com os slides: sim`.

Testado com Python 3.11, numpy 2.4 e SciPy 1.17.
