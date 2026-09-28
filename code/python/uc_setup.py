"""Acrescenta a sys.path as pastas do código da UC (Python).

Ao ser importado, acrescenta ao fim de sys.path, sem duplicar:
  * code/python/common (utilitários usados por vários decks);
  * todas as pastas sob code/python que contêm ficheiros .py (as pastas dos
    métodos), ignorando pastas escondidas e __pycache__.
É idempotente: importá-lo (ou chamar add_paths()) várias vezes não duplica
entradas.

Os exemplos ex*.py importam-no logo no início. Para usar as funções da UC nos
seus próprios scripts:
    import sys
    sys.path.insert(0, "<repo>/code/python")
    import uc_setup

Otimização — código da UC.  J. F. A. Madeira — Licença MIT.
"""
import os
import sys

RAIZ = os.path.dirname(os.path.abspath(__file__))


def _pastas():
    """common/ primeiro; depois as pastas com .py, por ordem alfabética."""
    comum = os.path.join(RAIZ, "common")
    pastas = [comum] if os.path.isdir(comum) else []
    for dirpath, dirnames, filenames in os.walk(RAIZ):
        dirnames[:] = sorted(d for d in dirnames
                             if not d.startswith(".") and d != "__pycache__")
        if dirpath != RAIZ and dirpath != comum and any(f.endswith(".py") for f in filenames):
            pastas.append(dirpath)
    return pastas


def add_paths():
    """Acrescenta as pastas a sys.path (sem duplicar); devolve a lista."""
    ja = {os.path.normcase(os.path.realpath(p)) for p in sys.path if p}
    pastas = _pastas()
    for p in pastas:
        chave = os.path.normcase(os.path.realpath(p))
        if chave not in ja:
            sys.path.append(p)
            ja.add(chave)
    return pastas


add_paths()
