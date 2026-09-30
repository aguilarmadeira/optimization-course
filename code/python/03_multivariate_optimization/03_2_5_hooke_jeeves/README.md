# 3.2.5 — Hooke–Jeeves (Python)

Estado: **disponível**.

- Módulos/funções: `exploratory.py` — `exploratory(f, xb, fb, P, menor=None)` → `ExploratoryResult`; `hooke_jeeves.py` — `hooke_jeeves(f, x0, a, P0, T, kmax=10000, verbose=False, menor=None)` → `HookeJeevesResult` (usa `exploratory`). A exploração avalia x_j + P_j **e** x_j − P_j e fica com o melhor dos três pontos, como em Deb (2012, sec. 3.3.3). `menor` (opcional, por omissão u < v) é a comparação usada; o deck 4.3.3 usa-a para a regra de admissibilidade.
- Exemplo: `ex03_2_5_hooke_jeeves.py`, que reproduz os números do deck 3.2.5 (os mesmos números, e a mesma saída, que o exemplo MATLAB).
- Apontamentos: `notes/pt/` (deck 3.2.5).

## Como correr

Só precisa de `numpy`. Na pasta desta secção:

```bash
python ex03_2_5_hooke_jeeves.py
```

## O que o exemplo imprime

A exploração à mão em f(x) = 3x₁² + x₂² − 12x₁ − 8x₂ a partir de (1,1); o exemplo completo (uma linha por movimento); o Rosenbrock a partir de (−1.5, 2) e o custo de cada movimento. No fim, a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Exploração à mão: x₁: (1.5, 1) → −18.25 e (0.5, 1) → −12.25, fica (1.5, 1); x₂: (1.5, 1.5) → −21 e (1.5, 0.5) → −15, fica (1.5, 1.5); 4 avaliações.
- Exemplo completo (a = 2, P₀ = (0.5, 0.5), T = (0.1, 0.1)): padrões (2,2) → (2, 2.5) aceite, (2.5, 3.5) → (2,4) aceite, (2, 5.5) → (2,5) rejeitado; a exploração em (2,4) falha com P = 0.5, 0.25, 0.125; para em (2,4), f = −28 (7 movimentos, n_f = 32).
- Rosenbrock (P₀ = 0.5, a = 2, T = 10⁻⁶, como no caderno da comparação): primeira avaliação com f < 10⁻⁴ em n_f = 430; no fim, 137 movimentos e n_f = 662.
- Custo: 4 (= 2n) avaliações por exploração; 5 (= 1 + 2n) por movimento de padrão.
- Última linha: `confere com os slides: sim`.

**Contagens.** `nfev` = 1 (f(x₀)) + 2n por exploração + 1 por ponto tentativo; f(x_b) nunca é reavaliado. Os resultados têm os campos `x`, `fx`, `nit`, `nfev`, `ngev`, `nhev`, `history`, `flag` (e `message`, `cols`, `ftrace`; `hooke_jeeves` também `P`), com as mesmas contagens que a versão MATLAB. `nit` é o número de movimentos (explorações a partir da base + movimentos de padrão).
Testado com Python 3.11 e numpy 2.4.
