# 4.3.3 — Restrições na prática: penalizar, rejeitar ou comparar (Python)

Estado: **disponível**.

> **Notação das aulas.** Nos decks 4.2 e 4.3 o multiplicador de uma desigualdade g_j ≥ 0 chama-se **λ_j** (e o de uma igualdade **β_ℓ**). No código mantém-se o nome `u` (`info.u`) para as estimativas de λ_j; os valores são os mesmos.

- Módulo/funções: `extreme_barrier.py` — `extreme_barrier(f, g, x)`: f(x) se todas as g_j(x) ≥ 0, senão `inf` sem avaliar f; `extreme_barrier(F, None, x)`: restrição oculta, `inf` se F falha (NaN, inf ou exceção).
- Usa o `hooke_jeeves` do deck 3.2.5 (com o argumento opcional `menor`, a comparação usada pela regra de admissibilidade) e o `nelder_mead` do 3.2.3, encontrados através de `import uc_setup`.
- Exemplo: `ex04_3_3_constraint_handling.py`, que reproduz os números do deck 4.3.3. Os mesmos números que o exemplo MATLAB.
- Apontamentos: `notes/pt/` e `notes/en/` (deck 4.3.3).

## Como correr

Só precisa de `numpy`. Na pasta desta secção:

```bash
python ex04_3_3_constraint_handling.py
```

## O que o exemplo imprime

O problema: min (x₁ − 2)² + (x₂ − 1)² s.a. g = 2 − x₁ − x₂ ≥ 0; x* = (1.5, 0.5), u* = 1. Hooke–Jeeves com x₀ = (0, 0), a = 2, P₀ = 0.5, T = 10⁻⁶.

1. A tabela dos cinco tratamentos (A ignorar, B quadrática com R = 1, 10, 100, 1000, C exata com R = 0.5, 1, 2, D barreira extrema, E regra de admissibilidade).
2. B: v(x_R) = 1/(1 + 2R).
3. Restrição oculta: um «simulador» que devolve NaN, tratado com `extreme_barrier(simulador, None, x)`.
4. D e E a partir de (3, 3), inadmissível.
5. O bloqueio em (1, 1): as quatro tentativas ±P e_j.
6. A correção: o mesmo HJ nas variáveis y, x = Q y, Q = [1 1; 1 −1]/√2.
7. Informativo: o Nelder–Mead do 3.2.3 em C (R = 2) e D.

## Saída esperada (resumo)

- A (2, 1), n_F = 96; B (1.667, 0.667), (1.524, 0.524), (1.502, 0.502), (1.500, 0.500) com n_F = 253, 244, 311, 644; C (1.75, 0.75), x*, (1, 1); D e E (1, 1) com n_F = 91 (D: 45 chamadas fora de X, sem avaliar f).
- Restrição oculta: a mesma corrida que D. De (3, 3): D não sai de x₀; E chega a (1, 1) com n_F = 96.
- Direções rodadas: C (R = 2), D e E chegam a x* com n_F = 208.
- Última linha: `confere com os slides: sim`.

**Contagens.** n_F = chamadas à função (ou ao par (v, f), na regra E) que o HJ recebe, incluindo f(x₀); na barreira extrema contam também as chamadas fora de X, em que f não é avaliada.
Testado com Python 3.11 e numpy 2.4.
