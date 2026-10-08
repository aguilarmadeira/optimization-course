# 3.2.5 — Hooke–Jeeves (Python)

Estado: **disponível**.

- Módulos/funções: `exploratory.py` — `exploratory(f, xb, fb, P, menor=None)` → `ExploratoryResult`; `hooke_jeeves.py` — `hooke_jeeves(f, x0, a, P0, T, kmax=10000, verbose=False, menor=None, variante="deb", R=2.0)` → `HookeJeevesResult` (usa `exploratory`). A exploração avalia x_j + Δ_j **e** x_j − Δ_j e fica com o melhor dos três pontos, como em Deb (2012, sec. 3.3.3). `menor` (opcional, por omissão u < v) é a comparação usada; o deck 4.3.3 usa-a para a regra de admissibilidade.
- Exemplo: `ex03_2_5_hooke_jeeves.py`, que reproduz os números do deck 3.2.5 (os mesmos números, e a mesma saída, que o exemplo MATLAB).
- Apontamentos: `notes/pt/` (deck 3.2.5).

O algoritmo é o das aulas (Deb, 2012, sec. 3.3.3), variante por omissão `deb`: explora-se em torno de x_k; se melhorar, faz-se o padrão x^P_{k+1} = 2x_k − x_{k−1} e explora-se em torno de x^P_{k+1}; se o resultado não for melhor do que x_k (ou se a exploração em torno de x_k falhar), passo 3: se ‖Δ‖ < ε, pára; senão, fazer Δ = Δ/R e explorar de novo em torno de x_k. Aqui `T` é ε (um escalar). A variante `1961` (Hooke e Jeeves, 1961) volta a explorar com o mesmo Δ depois de um padrão falhado e pára quando todos os Δ_j < T_j; é a usada no deck 4.3.3.

## Como correr

Só precisa de `numpy`. Na pasta desta secção:

```bash
python ex03_2_5_hooke_jeeves.py
```

## O que o exemplo imprime

A exploração à mão em f(x) = 3x₁² + x₂² − 12x₁ − 8x₂ a partir de (1,1); o exemplo completo (uma linha por movimento); o exemplo de Himmelblau do Deb; o Rosenbrock a partir de (−1.5, 2) e o custo de cada movimento. No fim, a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Exploração à mão (Δ = (0.5, 0.5)): x₁: f⁺ = f(1.5, 1) = −18.25 e f⁻ = f(0.5, 1) = −12.25, fica (1.5, 1); x₂: f⁺ = f(1.5, 1.5) = −21 e f⁻ = f(1.5, 0.5) = −15, fica (1.5, 1.5); 4 avaliações.
- Exemplo completo (Δ = (0.5, 0.5), R = 2, ε = 0.2): padrões x^P = (2,2), f = −24 → (2, 2.5), −25.75, sucesso; (2.5, 3.5), −27 → (2,4), −28, sucesso; (2, 5.5), −25.75 → (2,5), −27, insucesso (passo 3); a exploração em (2,4) falha com Δ = 0.25 e 0.125; ‖Δ‖ = 0.177 < ε: pára em (2,4), f = −28 (6 movimentos, n_f = 28).
- Himmelblau a partir de (0,0), Δ = (0.5, 0.5), ε = 0.2 (Deb, 2012, exercício 3.3.3): x₁ = (0.5, 0.5), x₂ = (1.5, 1.5), x₃ = (3,2); o padrão (4.5, 2.5) → (4,2), f = 50, falha; pára em (3,2), f = 0.
- Rosenbrock (Δ = 0.5, R = 2, ε = 10⁻⁶): primeira avaliação com f < 10⁻⁴ em n_f = 627; no fim, 140 movimentos e n_f = 680. Com a versão de 1961 (`variante` = 1961): n_f = 430.
- Custo: 4 (= 2n) avaliações por exploração; 5 (= 1 + 2n) por movimento em padrão.
- Última linha: `confere com os slides: sim`.

**Contagens.** `nfev` = 1 (f(x₀)) + 2n por exploração + 1 por ponto de padrão x^P; f(x_b) nunca é reavaliado. Os resultados têm os campos `x`, `fx`, `nit`, `nfev`, `ngev`, `nhev`, `history`, `flag` (e `message`, `cols`, `ftrace`; `hooke_jeeves` também `P`), com as mesmas contagens que a versão MATLAB. `nit` é o número de movimentos (explorações a partir da base + movimentos de padrão).
Testado com Python 3.11 e numpy 2.4.
