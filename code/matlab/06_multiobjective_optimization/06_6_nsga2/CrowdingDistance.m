function d = CrowdingDistance(F, idx)
%CROWDINGDISTANCE  Distância de aglomeração (crowding distance) numa frente.
%
%   d = CrowdingDistance(F, idx)
%
%   Os extremos recebem d = Inf; para os outros, ordenando por cada
%   objetivo i,
%       d_p = soma_i (f_i(seguinte) - f_i(anterior)) / (f_i^max - f_i^min).
%   Frentes com 1 ou 2 pontos: d = Inf para todos.
%
%   Entradas
%     F        matriz N x m dos objetivos
%     idx      índices dos pontos da frente
%
%   Saídas
%     d        vetor coluna com d_p, pela ordem de idx
%              (maior d = região menos povoada = maior prioridade)
%
%   Otimização — deck 6.6 (NSGA-II), frame «Distância de aglomeração».
%   Usada por NSGA2 (ver ex06_6_nsga2.m).
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas.
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

idx = idx(:);
m = size(F, 2);
L = numel(idx);
d = zeros(L, 1);
if L <= 2
  d(:) = Inf;
  return
end
for i = 1:m
  [v, o] = sort(F(idx, i));          % ordenar a frente pelo objetivo i
  amp = v(end) - v(1);               % f_i^max - f_i^min na frente
  d(o(1)) = Inf;  d(o(end)) = Inf;   % extremos
  if amp > 0
    for k = 2:L - 1
      d(o(k)) = d(o(k)) + (v(k + 1) - v(k - 1))/amp;
    end
  end
end
end
