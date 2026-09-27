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
