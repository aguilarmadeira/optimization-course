# 3.2.5 — Hooke–Jeeves (MATLAB)

Estado: **disponível**.

- Funções: `Exploratory.m` — `[xe, fe, info] = Exploratory(f, xb, fb, P, menor)`; `HookeJeeves.m` — `[x, fx, info] = HookeJeeves(f, x0, a, P0, T, kmax, verbose, menor)` (usa `Exploratory`). A exploração avalia x_j + P_j **e** x_j − P_j e fica com o melhor dos três pontos, como em Deb (2012, sec. 3.3.3). `menor` (opcional, por omissão `@(u,v) u < v`) é a comparação usada; o deck 4.3.3 usa-a para a regra de admissibilidade.
- Exemplo: `ex03_2_5_hooke_jeeves.m`, que reproduz os números do deck 3.2.5.
- Apontamentos: `notes/pt/` (deck 3.2.5).

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave:

```matlab
ex03_2_5_hooke_jeeves
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex03_2_5_hooke_jeeves"`). A ajuda: `help HookeJeeves`, `help Exploratory`.

## O que o exemplo imprime

A exploração à mão em f(x) = 3x₁² + x₂² − 12x₁ − 8x₂ a partir de (1,1); o exemplo completo (uma linha por movimento: exploração ou padrão, base, ponto tentativo, ponto explorado, P, resultado, n_f acumulado); o Rosenbrock a partir de (−1.5, 2) e o custo de cada movimento. No fim, a verificação com os slides. Os slides usam vírgula decimal; o código usa o ponto.

## Saída esperada (resumo)

- Exploração à mão: x₁: (1.5, 1) → −18.25 e (0.5, 1) → −12.25, fica (1.5, 1); x₂: (1.5, 1.5) → −21 e (1.5, 0.5) → −15, fica (1.5, 1.5); 4 avaliações.
- Exemplo completo (a = 2, P₀ = (0.5, 0.5), T = (0.1, 0.1)): padrão (2,2), f = −24 → (2, 2.5), −25.75, aceite; (2.5, 3.5), −27 → (2,4), −28, aceite; (2, 5.5), −25.75 → (2,5), −27, rejeitado; a exploração em (2,4) falha com P = 0.5, 0.25, 0.125; P = 0.0625 < T: para em (2,4), f = −28 (7 movimentos, n_f = 32).
- Rosenbrock (P₀ = 0.5, a = 2, T = 10⁻⁶, como no caderno da comparação): primeira avaliação com f < 10⁻⁴ em n_f = 430; no fim, 137 movimentos e n_f = 662. Com T = 10⁻⁵ o n_f até f < 10⁻⁴ é o mesmo (430).
- Custo: 4 (= 2n) avaliações por exploração; 5 (= 1 + 2n) por movimento de padrão.
- Última linha: `confere com os slides: sim`.

**Contagens.** `info.nfev` = 1 (f(x₀)) + 2n por exploração + 1 por ponto tentativo x_t. f(x_b) fica guardado e nunca é reavaliado; f(x) final não é reavaliado. `Exploratory` recebe f(x_b) já conhecido e não o conta. `info.ngev = info.nhev = 0`.
`info.nit` é o número de movimentos (explorações a partir da base + movimentos de padrão); `info.history` tem uma linha por movimento, `[m tipo xb xt f(xt) xe f(xe) f_ref P sucesso nfev]`. `info.ftrace` guarda os valores de f pela ordem de avaliação: `find(info.ftrace < 1e-4, 1)` dá as avaliações até f < 10⁻⁴.
P só diminui (P ← P/2 depois de uma exploração falhada em torno da base); não é reposto em P₀ a cada sucesso.
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave. O `patternsearch` do MATLAB (Global Optimization Toolbox) não é usado.
