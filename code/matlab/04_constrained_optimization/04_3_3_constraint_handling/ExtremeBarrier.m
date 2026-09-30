function y = ExtremeBarrier(f, g, x)
%EXTREMEBARRIER  Barreira extrema: +Inf fora da região admissível, sem avaliar f.
%
%   y = ExtremeBarrier(f, g, x)    restrições explícitas g(x) >= 0 (vetor;
%                                  convenção do cap. 4): y = f(x) se todas
%                                  as g_j(x) >= 0; senão
%                                  y = Inf e f NÃO é avaliada.
%   y = ExtremeBarrier(F, [], x)   restrição oculta (não há g): y = F(x), ou
%                                  Inf se a avaliação falha (NaN, Inf ou erro).
%
%   Uso: F = @(x) ExtremeBarrier(f, g, x); depois, qualquer método que só
%   compare valores (p. ex. HookeJeeves, 3.2.5) minimiza F. É o que fazem o
%   DMS e o MultiGLODS: f = +Inf nos pontos inadmissíveis, que não são
%   avaliados.
%
%   Otimização — deck 4.3.3 (restrições na prática).
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if ~isempty(g)
  if any(g(x) < 0), y = Inf; else, y = f(x); end
else
  try
    y = f(x);
  catch
    y = Inf;                     % o simulador falhou com erro
  end
  if isnan(y), y = Inf; end      % ... ou devolveu NaN
end
end
