function [fronts, rank] = NonDominatedSort(F, cv)
%NONDOMINATEDSORT  Ordenação não dominada (fast non-dominated sort) e regra de Deb.
%
%   [fronts, rank] = NonDominatedSort(F)
%   [fronts, rank] = NonDominatedSort(F, cv)
%
%   Entradas
%     F        matriz N x m dos valores dos objetivos (uma linha por ponto)
%     cv       (opcional, []) vetor com a violação das restrições de cada
%              ponto (cv = 0: admissível). Com cv, aplica a regra de Deb:
%              admissível vence inadmissível; entre inadmissíveis, menor
%              violação; entre admissíveis, dominância normal. Os admissíveis
%              são ordenados em frentes; cada inadmissível forma uma frente
%              sozinho, depois, por violação crescente (empates pela ordem
%              dos índices).
%
%   Saídas
%     fronts   cell array {F_1, F_2, ...}, cada um vetor linha de índices
%     rank     vetor coluna N x 1 com o rank de cada ponto (1 = frente F_1)
%
%   Otimização — deck 6.6 (NSGA-II), frame «Ordenação não dominada: o rank»;
%   dominância: deck 6.2. Usada por NSGA2 (ver ex06_6_nsga2.m).
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas (p. ex. não trata valores de f que sejam NaN).
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 2, cv = []; end
N = size(F, 1);
if isempty(cv)
  fronts = fastsort(F);
else
  cv = cv(:);
  adm = find(cv <= 0)';
  inadm = find(cv > 0)';
  fronts = {};
  if ~isempty(adm)
    fr = fastsort(F(adm, :));
    for k = 1:numel(fr)
      fronts{end + 1} = adm(fr{k});                  %#ok<AGROW>
    end
  end
  [~, o] = sort(cv(inadm));                          % sort é estável
  for i = inadm(o)
    fronts{end + 1} = i;                             %#ok<AGROW>
  end
end
rank = zeros(N, 1);
for k = 1:numel(fronts)
  rank(fronts{k}) = k;
end
end

function fronts = fastsort(F)
% Ordenação de Deb et al. (2002): S{p} = quem p domina, n(p) = quantos dominam p.
N = size(F, 1);
S = cell(N, 1);
n = zeros(N, 1);
fronts = {[]};
for p = 1:N
  Fp = repmat(F(p, :), N, 1);
  domina = all(Fp <= F, 2) & any(Fp < F, 2);       % p domina q
  dominado = all(F <= Fp, 2) & any(F < Fp, 2);     % q domina p
  S{p} = find(domina)';
  n(p) = sum(dominado);
  if n(p) == 0
    fronts{1}(end + 1) = p;                          % ninguém domina p: F_1
  end
end
i = 1;
while ~isempty(fronts{i})                            % retira F_i e repete
  nxt = [];
  for p = fronts{i}
    for q = S{p}
      n(q) = n(q) - 1;
      if n(q) == 0
        nxt(end + 1) = q;                            %#ok<AGROW>
      end
    end
  end
  i = i + 1;
  fronts{i} = nxt;
end
fronts = fronts(1:end - 1);
end
