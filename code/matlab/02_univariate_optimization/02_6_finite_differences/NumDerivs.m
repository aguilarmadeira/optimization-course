function [df, ddf, info] = NumDerivs(f, x, h)
%NUMDERIVS  Derivadas numéricas f'(x) e f''(x) por diferenças centradas.
%
%   [df, ddf, info] = NumDerivs(f, x)
%   [df, ddf, info] = NumDerivs(f, x, h)
%
%   Com os mesmos três valores f(x-h), f(x), f(x+h) (3 avaliações de f):
%     df  = (f(x+h) - f(x-h)) / (2h)            erro O(h^2)
%     ddf = (f(x+h) - 2 f(x) + f(x-h)) / h^2    erro O(h^2)
%
%   Entradas
%     f   função (handle)
%     x   ponto
%     h   (opcional) passo; por omissão a regra da UC
%         h = 0.01|x| se |x| > 0.01, 1e-4 caso contrário (= max(0.01|x|, 1e-4))
%
%   Saídas
%     df, ddf   aproximações centradas de f'(x) e f''(x)
%     info      estrutura com
%       .nfev     3 (f(x-h), f(x), f(x+h))
%       .ngev, .nhev   0 (derivadas numéricas contam como avaliações de f)
%       .h        passo usado
%       .fm, .f0, .fp   f(x-h), f(x), f(x+h)
%       .df_fwd   diferença progressiva (f(x+h) - f(x))/h, erro O(h)
%       .df_bwd   diferença regressiva  (f(x) - f(x-h))/h, erro O(h)
%                 (as duas reutilizam os mesmos valores: sem custo extra)
%
%   h grande: erro de truncatura; h pequeno: erro de arredondamento
%   (cancelamento). A regra de 1 % é uma regra de trabalho, não um passo
%   ótimo universal.
%
%   Otimização — deck 2.6 (Derivadas numéricas).
%   Reproduz o exemplo dos slides: ver ex02_6_finite_differences.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas.
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 3 || isempty(h)
  h = max(0.01*abs(x), 1e-4);
end

% --- núcleo (igual ao dos slides) -----------------------------------------
fm = f(x - h);  f0 = f(x);  fp = f(x + h);
df  = (fp - fm) / (2*h);
ddf = (fp - 2*f0 + fm) / h^2;
% -------------------------------------------------------------------------

info.nfev = 3;
info.ngev = 0;
info.nhev = 0;
info.h = h;
info.fm = fm;  info.f0 = f0;  info.fp = fp;
info.df_fwd = (fp - f0) / h;
info.df_bwd = (f0 - fm) / h;
end
