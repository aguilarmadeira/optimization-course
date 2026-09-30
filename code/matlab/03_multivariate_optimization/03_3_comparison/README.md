# 3.3 / 3.3.3 — Comparação dos métodos do capítulo 3 (MATLAB)

Estado: **disponível**.

- Função: `fg.m` — `[f, g, H] = fg(x)`: Rosenbrock, com o gradiente só quando `nargout > 1` e a Hessiana só quando `nargout > 2` (a função «completa» que o deck 3.3 refere para o `fminunc`).
- Contador: `contador.m` — conta n_f, n_g, n_H de uma corrida e guarda (n_f, n_g, n_H) na primeira avaliação com f < alvo.
- Exemplo: `ex03_3_comparison.m`, que reproduz a tabela «Os métodos estudados no mesmo problema, com o custo em avaliações» do deck 3.3.3.
- Usa as funções da UC das outras pastas do capítulo e `LineSearch` (de `common/`), postas no caminho por `uc_setup` (ver `code/README.md`): `LocalRandomSearch` (3.2.1), `NelderMead` (3.2.3), `BoxEvo` (3.2.4), `HookeJeeves` (3.2.5), `SteepestDescent` (3.3.1), `NewtonND` (3.3.2), `ConjGrad` (3.3.3).
- Apontamentos: `notes/pt/` (decks 3.3 e 3.3.3). Versão Python e caderno `cap3_comparacao.ipynb`: `code/python/03_multivariate_optimization/03_3_comparison/`.

## Como correr

Na pasta desta secção, em MATLAB ou GNU Octave (cerca de 20 s em Octave):

```matlab
ex03_3_comparison
```

(em Octave, na linha de comandos: `octave --no-gui --quiet --eval "ex03_3_comparison"`).

O `fminunc` do slide «Os contadores já feitos»:

```matlab
opts = optimoptions('fminunc', 'SpecifyObjectiveGradient', true);   % Octave: optimset('GradObj', 'on')
[x, fx, ~, out] = fminunc(@fg, [-1.5; 2], opts);
out.funcCount, out.iterations
```

## Regras da comparação

Mesmo problema (Rosenbrock), mesmo ponto inicial (−1.5, 2), mesmo critério: custo até à **primeira avaliação de f com f < 10⁻⁴** (f nos iterandos dos métodos com derivadas também conta em n_f). Custo em avaliações, n_f, n_g, n_H à parte; equivalente em f com derivadas centrais (n = 2): n_f + 4n_g + 4n_H. Métodos de 3.3 com a pesquisa em linha precisa (Brent, tol 10⁻¹⁰, `LineSearch`). Parâmetros: Nelder–Mead com o simplex x₀, x₀ + 0.5e₁, x₀ + 0.5e₂ (ε_x = 10⁻⁸, ε_f = 10⁻¹⁰); Box Δ₀ = (1,1), ε_x = 10⁻⁶; Hooke–Jeeves a = 2, P₀ = 0.5, T = 10⁻⁶; gradiente e FR com ε_g = 10⁻⁸, k_max = 3000 (FR com reinício a cada n); Newton ε_g = 10⁻⁸ (amortecido com `damped` e `modify`, como no caderno); pesquisa aleatória localizada r₀ = 1, m = 20, γ = 0.9, 30 corridas com `rng(0)`…`rng(29)`.

## Saída esperada (resumo)

| Método | n_f | n_g | n_H | equiv. em f | slides |
|---|---:|---:|---:|---:|---|
| Aleatória localizada (30 corridas) | mediana 1030 [850; 1217] | — | — | 1030 | 975 [867; 1052] — ver nota |
| Nelder–Mead | 166 | — | — | 166 | = |
| Box (Δ₀ = 1) | 13 445 | — | — | 13 445 | = |
| Hooke–Jeeves (P₀ = 0.5) | 430 | — | — | 430 | = |
| Gradiente (não atinge f < 10⁻⁴ em 3000 it.; f = 2.4·10⁻⁴) | 60 900 | 3000 | — | 72 900 | = |
| Grad. conjugados (FR, reinício n) | 718 | 34 | — | 854 | = |
| Newton puro | 6 | 5 | 5 | 46 | = |
| Newton amortecido | 280 | 12 | 12 | 376 | 281 / 377 — ver nota |

Também: até f < 10⁻⁴, FR tem n_f = 684 nas pesquisas em linha + 34 iterandos = 718 (o texto do deck 3.3.3). Última linha: `confere com os slides: sim`.

**Pesquisa aleatória localizada (estocástica).** O gerador do MATLAB/Octave não é o `numpy.random.default_rng` usado nos slides; por isso a verificação aqui é **estatística**: mediana a menos de 15 % de 975 e f < 10⁻⁴ atingido nas 30 corridas (em Octave 8.4: mediana 1030, quartis [850; 1217], 30/30). A versão Python reproduz exatamente 975 [867; 1052].

**Newton amortecido.** As iterações, n_g e n_H coincidem com os slides; o n_f até f < 10⁻⁴ depende do último bit da solução de H d = −g, porque a pesquisa de Brent com tol 10⁻¹⁰ perto do ótimo muda uma ou outra avaliação: o `numpy` (LU) dá 281, como os slides e a versão Python; o `H \ g` do Octave (Cholesky) dá 280 (com `lu` explícito 284, com a regra de Cramer 293). O exemplo aceita n_f a menos de 2 % e di-lo na saída.

**Contagens.** O `contador` conta as avaliações feitas pelos métodos; no fim de cada corrida o exemplo confirma que os contadores coincidem com `info.nfev`, `info.ngev`, `info.nhev` da própria função (se não coincidirem imprime um AVISO). A parte do `fminunc` (quasi-Newton, critério e pesquisa em linha do solver) é informativa e não entra na verificação.
Testado em GNU Octave 8.4; só usa funções comuns a MATLAB e Octave (o `fminunc` só é usado se existir).
