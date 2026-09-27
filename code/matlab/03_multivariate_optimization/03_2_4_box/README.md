# 3.2.4 — Método de Box (EVOP simplificado) (MATLAB)

Estado: **disponível**.

- Função: `BoxEvo.m` — `[x, fx, info] = BoxEvo(f, x0, Delta, tolx, kmax, verbose)`
- Exemplo: `ex03_2_4_box.m`, que reproduz os números do deck 3.2.4.
- Apontamentos: `notes/pt/` (deck 3.2.4).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex03_2_4_box
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex03_2_4_box"`). A ajuda da função: `help BoxEvo`.

## O que o exemplo imprime

A tabela das 20 iterações no Himmelblau a partir de (0,0) com Δ = (2,2) e ε_x = 10⁻⁴ (centro, Δ, f, decisão mover/encolher, n_f acumulado), os valores de f nos 4 vértices da «iteração à mão», o Rosenbrock a partir de (−1.5, 2) com Δ₀ = (1,1) e o teste das diagonais em (0,0). No fim, a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Iteração à mão: f(0,0) = 170; vértices (−1,−1), (1,−1), (−1,1), (1,1): 170, 146, 130, 106.
- Himmelblau: tabela k = 0..9 igual à do slide; 20 iterações, n_f = 4 × 20 + 1 = 81; termina em (3,2), f = 0.
- Rosenbrock (ε_x = 10⁻⁶, como no caderno da comparação; o deck não indica ε_x): primeira avaliação com f < 10⁻⁴ em n_f = 13 445, com Δ = 6.1·10⁻⁵ (≈ 6·10⁻⁵); no fim, 6474 iterações e n_f = 25 897.
- Diagonais em (0,0) com passo 0.25: f = 4.08 a 11.33 (nenhuma melhora f = 1); f(0.25, 0) = 0.953.
- Última linha: `confere com os slides: sim`.

**Contagens.** `info.nfev = 1 + 2^n · info.nit`: f(x₀) conta, e em cada iteração avaliam-se **sempre** os 2ⁿ vértices, sem reaproveitar pontos já avaliados (é a convenção das contagens do deck). f(x) final não é reavaliado. `info.ngev = info.nhev = 0`.
`info.history` tem uma linha por iteração k = 0, …, nit−1: `[k x_k Δ f(x_k) decisão f_novo nfev]` (decisão 1 = mover, 0 = encolher), como a tabela do slide. `info.ftrace` guarda os valores de f pela ordem de avaliação: `find(info.ftrace < 1e-4, 1)` dá as avaliações até f < 10⁻⁴.
Ordem dos vértices: a de `dec2bin(0:2^n-1)`, como no slide (em 2D: (−1,−1), (−1,1), (1,−1), (1,1)); em caso de empate, `min` escolhe o primeiro.
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave.
