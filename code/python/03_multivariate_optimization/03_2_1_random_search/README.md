# 3.2.1 — Pesquisa aleatória pura e localizada (Python)

Estado: **disponível**.

- Módulos/funções: `random_search.py` — `random_search(f, a, b, N, s=None, ftarget=-inf, verbose=False)` → `RandomSearchResult`; `local_random_search.py` — `local_random_search(f, x0, r0, m=20, gamma=0.9, tolx=1e-8, kmax=5000, ftarget=-inf, s=None, verbose=False)` → `LocalRandomSearchResult`
- Exemplo: `ex03_2_1_random_search.py`, que reproduz exatamente os números do deck 3.2.1.
- Apontamentos: `notes/pt/` (deck 3.2.1).

## Como correr

Só precisa de `numpy`. Na pasta desta secção:

```bash
python ex03_2_1_random_search.py
```

## Método estocástico: reprodução exata

`s` é a semente (um inteiro) ou um `numpy.random.Generator` já criado. As funções usam `numpy.random.default_rng(s)` e consomem os números aleatórios pela mesma ordem dos scripts das figuras do deck e do caderno do capítulo 3: `a + (b - a)*rng.random(n)` por amostra na pura, `x + r*(2*rng.random((m, n)) - 1)` por iteração na localizada. Com as mesmas sementes, os resultados são **os mesmos, número a número**:

- pura, ‖x‖² em [−5,5]²: semente 1 (uma corrida); gerador com semente 7 partilhado pelas 30 corridas de cada N (tabela da análise experimental);
- localizada no Rosenbrock: semente 3 (uma corrida); sementes 0–29 (30 corridas);
- pura vs. localizada com orçamento 3000: sementes 0–29.

O exemplo MATLAB/Octave não tem este gerador e faz uma verificação estatística (ver o README MATLAB).

## O que o exemplo imprime

1. ‖x‖² em [−5,5]², uma corrida, N = 1000: as melhorias do melhor ponto (`history`).
2. A tabela da análise experimental (N = 10, 100, 1000, 10 000; 30 corridas).
3. A localizada no Rosenbrock de (−1.5, 2), m = 20, r₀ = 1, γ = 0.9: uma corrida de 150 iterações e o custo até f < 10⁻⁴ em 30 corridas.
4. Pura vs. localizada com o mesmo orçamento (3000), 30 corridas.

## Saída esperada (resumo)

- N = 100: f = 1.406; N = 1000: f = 0.053 (n_f = N).
- Tabela: melhor f (mediana) 2.3 / 0.29 / 0.018 / 0.0028; ‖x_best‖ 1.5 [1.2; 2.1] / 0.54 [0.38; 0.76] / 0.14 [0.09; 0.17] / 0.053 [0.027; 0.070].
- Localizada, semente 3: 150 iterações, n_f = 3001 (com f(x₀)), x = (1, 1), f = 5.9·10⁻¹³.
- Custo até f < 10⁻⁴ (`nhit`), 30 corridas: 975 [867; 1052] (mediana [Q1; Q3]).
- Orçamento 3000: pura 0.014 [0.008; 0.022]; localizada 2·10⁻¹² [5·10⁻¹³; 8·10⁻¹²].
- Última linha: `confere com os slides: sim`.

O resultado tem os campos `x`, `fx`, `nit`, `nfev`, `nhit`, `ngev`, `nhev`, `history`, `flag` (e `message`), com as mesmas contagens que a versão MATLAB: na localizada `nfev = 1 + m·nit` (inclui f(x₀)); `nhit` é a primeira avaliação com f < `ftarget` (f(x₀) é a 1.ª; `inf` se nunca).
Testado com Python 3.11 e numpy 2.4.
