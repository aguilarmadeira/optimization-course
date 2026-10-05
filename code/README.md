# Código da UC — organização e caminhos

- **Clonar o repositório inteiro.** As pastas dependem umas das outras (p. ex. o exemplo de 3.2.2 usa o Nelder–Mead de 3.2.3), por isso uma pasta copiada sozinha não corre.
- **Os exemplos correm a partir de qualquer pasta.** Cada `ex*.m` começa por

  ```matlab
  addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC
  ```

  e cada `ex*.py` por

  ```python
  import os, sys
  sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
  import uc_setup  # noqa: F401,E402  (acrescenta as pastas do código da UC ao caminho)
  ```

- **Nos seus próprios scripts**, para usar as funções da UC:
  - MATLAB/Octave: `addpath('<repo>/code/matlab'); uc_setup`
  - Python: `import sys; sys.path.insert(0, '<repo>/code/python'); import uc_setup`

  `uc_setup` põe no caminho `common/` e todas as pastas dos métodos (sem pastas escondidas); pode chamar-se mais de uma vez.
- **Utilitários comuns** (usados por vários decks): `matlab/common/` e `python/common/`.
  - `LineSearch.m` / `line_search.py` — pesquisa em linha dos métodos de 3.3 (gradiente, Newton, gradientes conjugados e a comparação do cap. 3).
  - `NelderMeadComp.m` / `nelder_mead_comp.py` — a variante do Nelder–Mead usada na comparação de 4.3 (contração exterior aceite só com `<`); minimizador interno de 4.3.1 e 4.3.2.
- **Cada função existe uma só vez.** A função de um método fica na pasta do deck que a introduz (p. ex. `NelderMead` em `03_2_3_nelder_mead/`); quem a reutiliza noutro deck usa-a de lá, sem cópia. O que é usado por vários decks e não pertence a um só fica em `common/`.
- O caderno `python/03_multivariate_optimization/03_3_comparison/cap3_comparacao.ipynb` é autónomo (tem as suas próprias implementações).

## Correr no Colab ou no MATLAB Online (sem instalar nada)

Há um caderno por exemplo (`ex*.ipynb`, ao lado do `ex*.py`) e um caderno por capítulo em `python/notebooks/`, com todos os exemplos desse capítulo. A primeira célula clona o repositório no Colab (ou, se o caderno for aberto dentro do repositório, usa `code/python`); cada secção corre um exemplo e confere-o com os slides.

| Capítulo | Colab | MATLAB Online |
|---|---|---|
| 2 Otimização unidimensional | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aguilarmadeira/optimization-course/blob/main/code/python/notebooks/02_univariate_optimization.ipynb) | [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=aguilarmadeira/optimization-course&file=code/matlab/02_univariate_optimization/ex02_capitulo.m) |
| 3 Otimização multidimensional sem restrições | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aguilarmadeira/optimization-course/blob/main/code/python/notebooks/03_multivariate_optimization.ipynb) | [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=aguilarmadeira/optimization-course&file=code/matlab/03_multivariate_optimization/ex03_capitulo.m) |
| 4 Otimização com restrições | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aguilarmadeira/optimization-course/blob/main/code/python/notebooks/04_constrained_optimization.ipynb) | [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=aguilarmadeira/optimization-course&file=code/matlab/04_constrained_optimization/ex04_capitulo.m) |
| 5 Otimização global | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aguilarmadeira/optimization-course/blob/main/code/python/notebooks/05_global_optimization.ipynb) | [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=aguilarmadeira/optimization-course&file=code/matlab/05_global_optimization/ex05_capitulo.m) |
| 6 Otimização multiobjetivo | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aguilarmadeira/optimization-course/blob/main/code/python/notebooks/06_multiobjective_optimization.ipynb) | [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=aguilarmadeira/optimization-course&file=code/matlab/06_multiobjective_optimization/ex06_capitulo.m) |

Os `ex*.py` são a fonte: os cadernos são gerados por `python/notebooks/make_notebooks.py` (`python make_notebooks.py --executar` gera os dois tipos e guarda as saídas). Depois de alterar um exemplo, voltar a gerá-los.

## Correr no MATLAB Online

O botão **MATLAB Online** de cada capítulo (na tabela acima e no topo de cada capítulo do [README](../README.pt.md)) abre o script `matlab/<capítulo>/ex0N_capitulo.m`, que corre, por ordem, todos os exemplos MATLAB desse capítulo.

- **É preciso conta MathWorks.**
- O link **clona o repositório para o MATLAB Drive** e abre o script do capítulo.
- Carrega-se em **Run**. Se o MATLAB perguntar, escolhe-se **«Change Folder»**. Cada exemplo começa com uma linha de separação e termina com «confere com os slides: sim» (exceto `ex06_6_gamultiobj`, uma variante com o `gamultiobj` que não tem verificação).
- **Toolboxes.** Alguns exemplos usam funções de toolboxes: `fmincon`/`linprog` (Optimization Toolbox) em `ex04_1_lagrange`, `ex04_2_kkt`, `ex06_3_aggregation` e `ex06_4_epsilon_constraint`, e `gamultiobj` (Global Optimization Toolbox) em `ex06_6_gamultiobj`. O `fminunc` de `ex03_3_comparison` está numa parte informativa protegida por `try`/`catch`. No MATLAB Online estes exemplos só correm se a licença incluir a toolbox; se não incluir, o script do capítulo diz qual exemplo não correu e continua com os seguintes.
- **Sem conta MathWorks**, descarrega-se a pasta `code/matlab/` (ou o repositório inteiro) e corre-se no MATLAB ou no GNU Octave: `ex0N_capitulo` corre o capítulo todo e cada `ex*.m` corre um exemplo. No GNU Octave, os exemplos com `fmincon`/`linprog` usam o `sqp` e o `glpk`, e `ex06_6_gamultiobj` só deixa um aviso (o Octave não tem `gamultiobj`).
