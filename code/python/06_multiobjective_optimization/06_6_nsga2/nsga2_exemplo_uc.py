"""NSGA-II (pymoo) no exemplo convexo da UC (capítulo 6).

    min f1(x) = x1^2 + x2^2
    min f2(x) = (x1 - 1)^2 + x2^2
    s.a. x1^2 + x2^2 <= 1,  0 <= x1, x2 <= 1

Instalação:  pip install pymoo matplotlib
Execução:    python nsga2_exemplo_uc.py
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from pymoo.core.problem import ElementwiseProblem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.indicators.hv import HV


class ExemploUC(ElementwiseProblem):
    def __init__(self):
        super().__init__(n_var=2, n_obj=2, n_ieq_constr=1, xl=[0.0, 0.0], xu=[1.0, 1.0])

    def _evaluate(self, x, out, *args, **kwargs):
        f1 = x[0] ** 2 + x[1] ** 2
        f2 = (x[0] - 1) ** 2 + x[1] ** 2
        out["F"] = [f1, f2]
        out["G"] = [x[0] ** 2 + x[1] ** 2 - 1]      # restrição na forma g(x) <= 0


problem = ExemploUC()
algorithm = NSGA2(pop_size=40)                     # N = 40 (SBX + mutação polinomial por omissão)

res = minimize(problem, algorithm, ("n_gen", 50),  # G = 50 gerações
               seed=1, verbose=False)              # semente fixa: resultado reprodutível

X, F = res.X, res.F                                # soluções (decisões) e objetivos da aproximação final
hv = HV(ref_point=np.array([1.1, 1.65]))(F)        # mesmo ponto de referência dos slides
hv_txt = f"{hv:.3f}".replace(".", ",")                 # vírgula decimal (PT)
print(f"{len(F)} pontos não dominados; avaliações: {res.algorithm.evaluator.n_eval}; HV = {hv_txt} (frente exata: 1,648)")

t = np.linspace(0, 1, 200)
plt.plot(t ** 2, (1 - t) ** 2, "-", color="0.7", lw=2, label="frente exata")
plt.plot(F[:, 0], F[:, 1], "o", color="#BE1E2D", ms=5, label="NSGA-II (pymoo)")
plt.xlabel("$f_1$"); plt.ylabel("$f_2$"); plt.legend(frameon=False); plt.gca().set_aspect("equal")
virgula = FuncFormatter(lambda v, pos: f"{v:g}".replace(".", ","))   # vírgula decimal nos eixos
plt.gca().xaxis.set_major_formatter(virgula); plt.gca().yaxis.set_major_formatter(virgula)
plt.title(f"NSGA-II, N = 40, 50 gerações, seed = 1: HV = {hv_txt}")
plt.savefig("nsga2_exemplo_uc.png", dpi=150, bbox_inches="tight")
plt.show()
