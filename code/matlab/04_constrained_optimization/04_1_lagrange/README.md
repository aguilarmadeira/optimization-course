# 4.1 — Multiplicadores de Lagrange (MATLAB)

Estado: **disponível**.

- Não há função da UC neste deck: o código é um exemplo.
- Exemplo: `ex04_1_lagrange.m`, que reproduz os números do deck 4.1:
  lata: r* = 3,745, h* = 7,49, A* ≈ 264,4 cm², λ* = −0,534 (dA*/dV = −λ* = 0,534).
- Solver: `fmincon` (`lambda.eqnonlin`, Optimization Toolbox), só se existir; em GNU Octave, `sqp` (Octave base).
- Apontamentos: `notes/pt/` (deck 4.1). O deck 4.1 remete para o 4.2, «No computador».

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex04_1_lagrange
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex04_1_lagrange"`).

## O que o exemplo faz

**(a) Pelas condições de Lagrange**, com as contas do slide, e verifica numericamente ∇f(x*) + Σ λ_ℓ ∇h_ℓ(x*) = 0 e h_ℓ(x*) = 0, na convenção das aulas 𝓛 = f + Σ λ_ℓ h_ℓ (com 𝓛 = f − Σ λ_ℓ h_ℓ, como em Nocedal e Wright, só muda o sinal de λ):

1. O exemplo resolvido, min 2x₁² + x₂² s.a. x₁ + x₂ = 1 (sistema linear em (x, λ)); inclui «Geometria (1)» (∇f(0,8; 0,2) = (3,2; 0,4), ∇fᵀd = 2,8) e «O que significa λ?» (declive de f*(δ) em δ = 0 por diferenças centrais = −λ*); e a Hessiana orlada, det = −6 < 0 (mínimo).
2. Candidatos, min x₁ + x₂ s.a. x₁² + x₂² = 2: (−1,−1) com λ = ½ (mínimo) e (1,1) com λ = −½ (máximo), classificados pela 2.ª ordem em T (dᵀ∇²𝓛 d = 2λ) e pela Hessiana orlada (det = −8 e 8).
3. A lata (V = 330 cm³): λ = −2/r, h = 2r, r* = (V/2π)^{1/3}; e a 2.ª ordem com d = (1,−4) ∈ T: dᵀ∇²𝓛 d = 12π > 0 (não está no slide; vem da auditoria).
4. Várias restrições: x* = (3/2, 1/2, 1), λ* = (−2, −1), f* = 7/2.
5. Quando a condição de Lagrange falha: ∇h₁, ∇h₂ dependentes em (1,0,0); o melhor λ (mínimos quadrados) deixa resíduo 1 — não há multiplicadores.

**(b) Com um solver** (os problemas 1, 3 e 4):

- **MATLAB com a Optimization Toolbox** (`exist('fmincon','file')`): `fmincon` com `nonlcon = @(x) deal([], h(x))`. O `fmincon` escreve 𝓛 = f + λᵀc_eq, como nas aulas, logo **λ = `lambda.eqnonlin`**, sem troca de sinal.
- **GNU Octave** (não tem `fmincon`): `sqp(x0, f, h, [])`. O `sqp` usa h(x) = 0 e 𝓛 = f − λᵀh, logo **λ das aulas = −`lambda`** (verificado: no exemplo resolvido o `sqp` dá +4/3 e o script mostra −4/3).
- Os códigos de saída do `sqp` 101 (convergência normal) e 104 (passo demasiado pequeno) são ambos paragem normal nestes problemas.

## Saída esperada (resumo)

- Exemplo resolvido: x* = (0.3333, 0.6667), λ* = −1.3333, f* = 0.6667; declive = 1.333333 = −λ*; Hessiana orlada: det = −6.
- Candidatos: (−1,−1), λ = 0.5, f = −2, det orlada = −8, mínimo; (1,1), λ = −0.5, f = 2, det orlada = 8, máximo.
- Lata: r* = 3.7449 cm, h* = 7.4899 cm, A* = 264.3568 cm², λ* = −0.5341 cm²/cm³.
- Várias restrições: x* = (1.5, 0.5, 1), λ* = (−2, −1), f* = 3.5.
- Lagrange falha: característica 1, resíduo mínimo 1.
- Solver (Octave `sqp`): os mesmos valores (exemplo resolvido: 6 it., 9 aval.; lata: 10 it., 15 aval.; várias restrições: 3 it., 6 aval.).
- Última linha: `confere com os slides: sim`.

Testado em GNU Octave 8.4 (ramo `sqp`); o ramo `fmincon` usa só a sintaxe documentada (`optimset`, `deal`) mas não foi corrido aqui (sem MATLAB). Só usa funções comuns a MATLAB e Octave.
