function Q = cruzarBLX(pais, pc, beta)
%CRUZARBLX  Cruzamento BLX-beta aos pares (auxiliar de GeneticAlgorithm).
%
%   Q = cruzarBLX(pais, pc)
%   Q = cruzarBLX(pais, pc, beta)
%
%   Os pais (linhas de pais, N x n) cruzam-se aos pares 1-2, 3-4, ..., cada
%   par com probabilidade pc:
%       x = a.*x_a + (1 - a).*x_b,   a ~ U(-beta, 1 + beta) por coordenada,
%   e o segundo descendente troca os papéis de x_a e x_b. Com beta > 0 os
%   descendentes podem sair do segmento entre os pais. Sem cruzamento, os
%   descendentes são cópias dos pais. Com N ímpar, o último pai passa sem par.
%   beta é opcional (0.5, BLX-0.5 como no slide).
%
%   Otimização — deck 5.3 (Algoritmos genéticos).
%   Implementação didática.
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 3 || isempty(beta), beta = 0.5; end
[N, n] = size(pais);
Q = pais;
for j = 1:2:N-1
  if rand < pc
    a = -beta + (1 + 2*beta)*rand(1, n);       % a ~ U(-beta, 1+beta)
    Q(j, :)     = a.*pais(j, :) + (1 - a).*pais(j + 1, :);
    Q(j + 1, :) = (1 - a).*pais(j, :) + a.*pais(j + 1, :);
  end
end
end
