# 3.3 / 3.3.3 — Comparação dos métodos do capítulo 3 (Python)

Estado: **disponível**.

- Módulo/função: `fg.py` — `fg(x, nout=1)`: Rosenbrock; `nout = 2` devolve também o gradiente e `nout = 3` a Hessiana (o equivalente do `fg.m` com `nargout`).
- Exemplo: `ex03_3_comparison.py`, que reproduz a tabela «Os métodos estudados no mesmo problema, com o custo em avaliações» do deck 3.3.3 com as funções da UC; os mesmos números que o exemplo MATLAB (ver as notas).
- Caderno: `cap3_comparacao.ipynb` — o caderno citado nos decks 3.3 e 3.3.3; faz a mesma comparação com implementações próprias dentro do caderno (as dos scripts das figuras) e mostra a figura custo × f e a parte com `scipy.optimize.minimize`.
- Usa as funções da UC das outras pastas do capítulo e `line_search` (de `common/`), acrescentadas ao `sys.path` por `uc_setup` (ver `code/README.md`): `local_random_search` (3.2.1), `nelder_mead` (3.2.3), `box_evo` (3.2.4), `hooke_jeeves` (3.2.5), `steepest_descent` (3.3.1), `newton_nd` (3.3.2), `conj_grad` (3.3.3).
- Apontamentos: `notes/pt/` (decks 3.3 e 3.3.3).

## Como correr

O exemplo só precisa de `numpy` (a parte final, informativa, usa o SciPy se estiver instalado). Na pasta desta secção:

```bash
python ex03_3_comparison.py
```

O caderno abre-se com Jupyter (precisa de `numpy`, `scipy` e `matplotlib`) e dá a mesma tabela que o exemplo. O caderno é autónomo: tem as suas próprias implementações e só importa `numpy`, `scipy` e `matplotlib` (não usa `uc_setup` nem outras pastas do repositório).

## Regras da comparação

Mesmo problema (Rosenbrock), mesmo ponto inicial (−1.5, 2), mesmo critério: custo até à **primeira avaliação de f com f < 10⁻⁴** (f nas iteradas dos métodos com derivadas também conta em n_f). Custo em avaliações, n_f, n_g, n_H à parte; equivalente em f com derivadas centrais (n = 2): n_f + 4n_g + 4n_H. Métodos de 3.3 com a pesquisa em linha precisa (Brent, tol 10⁻¹⁰, `line_search`). Parâmetros: Nelder–Mead com o simplex x₀, x₀ + 0.5e₁, x₀ + 0.5e₂ (ε_x = 10⁻⁸, ε_f = 10⁻¹⁰); Box Δ₀ = (1,1), ε_x = 10⁻⁶; Hooke–Jeeves a = 2, P₀ = 0.5, T = 10⁻⁶; gradiente e FR com ε_g = 10⁻⁸, k_max = 3000 (FR com reinício a cada n); Newton ε_g = 10⁻⁸ (amortecido com `damped=True, modify=True`, como no caderno); pesquisa aleatória localizada r₀ = 1, m = 20, γ = 0.9, 30 corridas com as sementes 0…29 (`numpy.random.default_rng`, as dos scripts das figuras).

## Saída esperada (resumo)

| Método | n_f | n_g | n_H | equiv. em f |
|---|---:|---:|---:|---:|
| Aleatória localizada (30 corridas) | 975 [867; 1052] | — | — | 975 |
| Nelder–Mead¹ | 244 | — | — | 244 |
| Box (Δ₀ = 1) | 13 445 | — | — | 13 445 |
| Hooke–Jeeves¹ (Δ = 0.5) | 627 | — | — | 627 |
| Gradiente (não atinge f < 10⁻⁴ em 3000 it.; f = 2.4·10⁻⁴) | 60 900 | 3000 | — | 72 900 |
| Grad. conjugados (FR, reinício n) | 718 | 34 | — | 854 |
| Newton puro | 6 | 5 | 5 | 46 |
| Newton amortecido | 281 | 12 | 12 | 377 |

¹ Versões das aulas (Deb, 2012): simplex com análise do erro Q (`eps = 1e-6`) e Hooke–Jeeves com redução de Δ logo após um padrão falhado. Com as versões padrão (sem f(x_c), como o `fminsearch`; Hooke–Jeeves de 1961, `variante = "1961"`): 166 e 430.

Todos iguais aos slides (e ao caderno). Também: até f < 10⁻⁴, Newton amortecido n_f = 269 nas pesquisas em linha + 12 iteradas = 281; FR 684 + 34 = 718 (o texto do deck 3.3.3). Última linha: `confere com os slides: sim`.

**Sensibilidade do Newton amortecido.** O n_f até f < 10⁻⁴ depende do último bit da solução de H d = −g (a pesquisa de Brent com tol 10⁻¹⁰ perto do ótimo muda uma ou outra avaliação): `np.linalg.solve` (LU) dá 281, como os slides; o `H \ g` do Octave (Cholesky) dá 280. Por isso a versão MATLAB aceita n_f a menos de 2 %.
**Contagens.** O `Contador` do exemplo conta as avaliações feitas pelos métodos; no fim de cada corrida o exemplo confirma que coincidem com `nfev`, `ngev`, `nhev` da própria função (se não coincidirem imprime um AVISO). A parte com o SciPy (critério e pesquisa em linha do SciPy, números diferentes) é informativa, como no caderno.
Testado com Python 3.11 e numpy 2.4 (SciPy 1.17 na parte informativa).
