"""Gera os cadernos Jupyter/Colab (um por capítulo) a partir dos exemplos ex*.py.

Os ficheiros ex*.py são a única fonte: cada caderno tem uma célula de
preparação (no Colab clona o repositório; dentro do repositório usa a pasta
code/python) e, por método, uma célula de texto (deck, funções usadas e a
descrição do exemplo) seguida do código do exemplo, sem as 3 linhas de
caminhos (substituídas pela célula de preparação).

Depois de alterar um exemplo, correr outra vez (a partir de qualquer pasta):
    python make_notebooks.py            # gera os cadernos sem saídas
    python make_notebooks.py --executar # gera e executa (guarda as saídas)

Otimização — código da UC.  J. F. A. Madeira — Licença MIT.
"""
import glob
import os
import re
import sys

import nbformat
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook

AQUI = os.path.dirname(os.path.abspath(__file__))
PY = os.path.dirname(AQUI)                     # code/python
REPO = "aguilarmadeira/optimization-course"
BLOB = f"https://github.com/{REPO}/blob/main/code/python"
COLAB = f"https://colab.research.google.com/github/{REPO}/blob/main/code/python/notebooks"

CAPITULOS = {
    "02_univariate_optimization": "2 — Otimização unidimensional",
    "03_multivariate_optimization": "3 — Otimização multidimensional sem restrições",
    "04_constrained_optimization": "4 — Otimização com restrições",
    "05_global_optimization": "5 — Otimização global",
    "06_multiobjective_optimization": "6 — Otimização multiobjetivo",
}

# Títulos das secções = subtítulos dos decks (por nome do exemplo).
TITULOS = {
    "ex02_2_golden_section": "2.2 Método da secção áurea",
    "ex02_3_powell": "2.3 Interpolação quadrática sucessiva (método de Powell)",
    "ex02_4_bisection": "2.4 Método da bisseção",
    "ex02_5_newton": "2.5 Método de Newton–Raphson",
    "ex02_5_newton_safeguarded": "2.5 Método de Newton–Raphson salvaguardado",
    "ex02_6_finite_differences": "2.6 Aproximação numérica de derivadas",
    "ex03_2_1_random_search": "3.2.1 Método de pesquisa aleatória",
    "ex03_2_2_grid_search": "3.2.2 Pesquisa em grelha",
    "ex03_2_3_nelder_mead": "3.2.3 Método do simplex de Nelder–Mead",
    "ex03_2_4_box": "3.2.4 Método de Box — versão simplificada de EVOP",
    "ex03_2_5_hooke_jeeves": "3.2.5 Método de Hooke–Jeeves",
    "ex03_3_1_steepest_descent": "3.3.1 Método do gradiente",
    "ex03_3_2_newton": "3.3.2 Método de Newton",
    "ex03_3_3_conjugate_gradients": "3.3.3 Método dos gradientes conjugados (Fletcher–Reeves)",
    "ex03_3_comparison": "3.3.3 Os métodos do capítulo 3 no mesmo problema",
    "ex04_1_lagrange": "4.1 Extremos condicionados: multiplicadores de Lagrange",
    "ex04_2_kkt": "4.2 Condições de Karush–Kuhn–Tucker (KKT)",
    "ex04_3_1_exterior_penalty": "4.3.1 Método da função de penalização exterior",
    "ex04_3_2_barrier": "4.3.2 Método da função de barreira / penalização interior",
    "ex05_2_simulated_annealing": "5.2 Simulated annealing",
    "ex05_3_genetic_algorithms": "5.3 Algoritmos genéticos",
    "ex06_3_aggregation": "6.3 Métodos de agregação de objetivos",
    "ex06_4_epsilon_constraint": "6.4 Método do ε-constrangimento",
    "ex06_6_nsga2": "6.6 NSGA-II (algoritmos genéticos multiobjetivo)",
}

# Secções sem exemplo em Python (código de investigação distribuído na aula).
SEM_EXEMPLO = {
    "05_global_optimization": [
        ("5.1 GLODS — Global and Local Optimization using Direct Search",
         "05_global_optimization/05_1_glods/README.md")],
    "06_multiobjective_optimization": [
        ("6.5 Direct MultiSearch (DMS)",
         "06_multiobjective_optimization/06_5_dms/README.md")],
}

CABECALHO = re.compile(
    r"^import os, sys\n"
    r"sys\.path\.insert\(0, os\.path\.join\(os\.path\.dirname\(os\.path\.abspath\(__file__\)\), \"\.\.\", \"\.\.\"\)\)\n"
    r"import uc_setup[^\n]*\n", re.M)

PREPARACAO = '''\
# Preparação: põe as funções da UC no caminho do Python.
# No Colab (ou noutra pasta), clona o repositório; dentro do repositório, usa code/python.
import os, sys
if os.path.isfile(os.path.join("..", "uc_setup.py")):
    RAIZ = os.path.abspath("..")
else:
    if not os.path.isdir("optimization-course"):
        !git clone -q --depth 1 https://github.com/aguilarmadeira/optimization-course
    RAIZ = os.path.abspath(os.path.join("optimization-course", "code", "python"))
if RAIZ not in sys.path:
    sys.path.insert(0, RAIZ)
import uc_setup  # acrescenta as pastas dos métodos e common/
print("Código da UC pronto.")'''


def modulos_da_uc(fonte):
    """Ficheiros .py da UC importados pelo exemplo (caminho relativo a code/python)."""
    fich = []
    for mod in re.findall(r"^from (\w+) import", fonte, re.M):
        achados = [p for p in glob.glob(os.path.join(PY, "**", mod + ".py"), recursive=True)
                   if os.sep + "notebooks" + os.sep not in p]
        if achados:
            fich.append(os.path.relpath(achados[0], PY).replace(os.sep, "/"))
    return fich


def seccao(caminho):
    nome = os.path.splitext(os.path.basename(caminho))[0]
    rel = os.path.relpath(caminho, PY).replace(os.sep, "/")
    fonte = open(caminho, encoding="utf-8").read()
    m = re.match(r'"""(.*?)"""\n', fonte, re.S)
    doc, corpo = (m.group(1), fonte[m.end():]) if m else ("", fonte)
    corpo, n = CABECALHO.subn("", corpo, count=1)
    if n != 1:
        sys.exit(f"{rel}: cabeçalho de caminhos não encontrado")
    # Descrição: sem a linha «Correr (de qualquer pasta)…» nem a assinatura final.
    linhas = [l for l in doc.strip().splitlines()
              if not l.startswith("Correr (de qualquer pasta)")
              and not l.startswith("Otimização — ")]
    descricao = "\n".join(linhas).strip()
    funcs = ", ".join(f"[`{os.path.basename(f)}`]({BLOB}/{f})" for f in modulos_da_uc(corpo))
    texto = f"## {TITULOS[nome]}\n\nExemplo [`{os.path.basename(rel)}`]({BLOB}/{rel})"
    if funcs:
        texto += f" · funções da UC: {funcs}"
    texto += f"\n\n```text\n{descricao}\n```"
    return [new_markdown_cell(texto), new_code_cell(corpo.strip() + "\n")]


def caderno(pasta, titulo):
    num = titulo.split(" ")[0]
    ficheiro = f"{pasta}.ipynb"
    cel = [new_markdown_cell(
        f"# Otimização — capítulo {num}: {titulo.split(' — ', 1)[1]}\n"
        f"### Exemplos em Python\n\n"
        f"[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)]({COLAB}/{ficheiro})\n\n"
        "Cada secção corre o exemplo `ex….py` do método e compara os resultados com os slides "
        "(no fim: *confere com os slides: sim*). Comece pela célula **Preparação**; no Colab, "
        "*Runtime → Run all* corre tudo.\n\n"
        "Os slides usam vírgula decimal; aqui usa-se o ponto. Este caderno é gerado a partir dos "
        "ficheiros `ex….py` por `make_notebooks.py`.\n\n"
        f"*Examples for Chapter {num} of the course \"Optimization\" (ISEL, J. F. A. Madeira). "
        "Code comments and printed messages are in Portuguese; the English slides are in "
        "[`notes/en/`](https://github.com/aguilarmadeira/optimization-course/tree/main/notes/en).*"),
        new_markdown_cell("## Preparação"),
        new_code_cell(PREPARACAO)]
    exemplos = sorted(glob.glob(os.path.join(PY, pasta, "*", "ex*.py")))
    exemplos = [e for e in exemplos if os.path.splitext(os.path.basename(e))[0] in TITULOS]
    blocos = [(os.path.basename(os.path.dirname(e)), seccao(e)) for e in exemplos]
    for tit, readme in SEM_EXEMPLO.get(pasta, []):
        blocos.append((os.path.basename(os.path.dirname(readme)), [new_markdown_cell(
            f"## {tit}\n\nCódigo de investigação do docente, distribuído na aula; não faz parte "
            f"do código da UC nem deste caderno. Ver [`README.md`]({BLOB}/{readme}).")]))
    if pasta == "03_multivariate_optimization":
        nb3 = "03_multivariate_optimization/03_3_comparison/cap3_comparacao.ipynb"
        blocos.append(("03_3_comparison~", [new_markdown_cell(
            "## Caderno da comparação do capítulo 3\n\n"
            f"O caderno autónomo [`cap3_comparacao.ipynb`]({BLOB}/{nb3}) faz a mesma comparação com "
            "gráficos. [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)]"
            f"(https://colab.research.google.com/github/{REPO}/blob/main/code/python/{nb3})")]))
    for _, c in sorted(blocos, key=lambda b: b[0]):
        cel += c
    nb = new_notebook(cells=cel, metadata={
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
        "colab": {"provenance": [], "toc_visible": True}})
    return os.path.join(AQUI, ficheiro), nb


def main():
    executar = "--executar" in sys.argv
    for pasta, titulo in CAPITULOS.items():
        destino, nb = caderno(pasta, titulo)
        if executar:
            from nbconvert.preprocessors import ExecutePreprocessor
            ExecutePreprocessor(timeout=600, kernel_name="python3").preprocess(
                nb, {"metadata": {"path": AQUI}})
        nbformat.validate(nb)
        nbformat.write(nb, destino)
        print(f"{os.path.basename(destino)}: {sum(c.cell_type == 'code' for c in nb.cells) - 1} exemplos")


if __name__ == "__main__":
    main()
