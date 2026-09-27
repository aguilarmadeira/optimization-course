function [g, nfev] = GradFD(f, x, h)
%GRADFD  Gradiente numérico por diferenças centrais, coordenada a coordenada.
%
%   [g, nfev] = GradFD(f, x)
%   [g, nfev] = GradFD(f, x, h)
%
%     df/dx_i ~ (f(x + h_i e_i) - f(x - h_i e_i)) / (2 h_i)      erro O(h_i^2)
%
%   Entradas
%     f   função (handle) de R^n em R
%     x   ponto (vetor com n componentes)
%     h   (opcional) passo, escalar ou vetor com n componentes; por omissão
%         a regra de 2.6 em cada coordenada: h_i = max(0.01|x_i|, 1e-4)
%
%   Saídas
%     g     gradiente aproximado (com a forma de x)
%     nfev  avaliações de f: 2n (conta em n_f; n_g = 0)
%
%   Mesmo compromisso de 2.6: h grande, erro de truncatura; h pequeno, erro
%   de arredondamento.
%
%   Otimização — decks 3.3 (Derivadas numéricas em n dimensões) e 3.3.1.
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas.
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

n = numel(x);
if nargin < 3 || isempty(h)
  h = max(0.01*abs(x(:)), 1e-4);
end
if isscalar(h), h = h*ones(n, 1); end
g = zeros(size(x));
nfev = 0;
% --- núcleo (o dos slides) -------------------------------------------------
for i = 1:n
  e = zeros(size(x));  e(i) = h(i);
  g(i) = (f(x + e) - f(x - e)) / (2*h(i));
  nfev = nfev + 2;
end
% -------------------------------------------------------------------------
end
