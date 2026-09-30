function [xe, fe, info] = Exploratory(f, xb, fb, P, menor)
%EXPLORATORY  Movimento exploratório de Hooke–Jeeves em torno de xb.
%
%   [xe, fe, info] = Exploratory(f, xb, fb, P)
%   [xe, fe, info] = Exploratory(f, xb, fb, P, menor)
%
%   Percorre as coordenadas j = 1, ..., n, sequencialmente: avalia x_j + P_j e
%   x_j - P_j e fica com o melhor dos três pontos (o atual e as duas
%   tentativas), como em Deb (2012, sec. 3.3.3); se nenhuma tentativa melhora,
%   deixa x_j como estava. Cada coordenada parte do melhor ponto até aí.
%
%   Entradas
%     f     função (handle) de um vetor linha x (1 x n)
%     xb    ponto base (vetor; é tratado como linha)
%     fb    f(xb), já conhecido (NÃO é avaliado aqui nem contado)
%     P     perturbações, uma por coordenada (escalar = igual em todas)
%     menor (opcional) handle menor(u, v): o valor u é melhor do que v?
%           Por omissão, u < v. O método só compara valores, por isso pode
%           usar outra ordem: p. ex. a regra de admissibilidade (deck 4.3.3),
%           com f a devolver [v(x) f(x)]; .history e .ftrace registam o
%           último componente.
%
%   Saídas
%     xe    ponto explorado (vetor linha); xe = xb se a exploração falhou
%     fe    f(xe)
%     info  estrutura com
%       .nfev     avaliações de f: 2n (duas tentativas por coordenada)
%       .ngev, .nhev   0
%       .nit      coordenadas percorridas (= n)
%       .history  uma linha por avaliação: [j  sinal  xp(1:n)  f(xp)  aceite]
%                 sinal = +1 ou -1 (x_j + P_j ou x_j - P_j), aceite = 1 se foi
%                 a tentativa escolhida (empate entre +P_j e -P_j: fica +P_j)
%       .cols     nomes das colunas de .history
%       .ftrace   valores de f pela ordem em que foram avaliados (1 x nfev)
%       .flag     0 se melhorou (fe < fb); 1 se falhou (xe = xb)
%       .message  mensagem
%
%   Otimização — deck 3.2.5 (Método de Hooke–Jeeves).
%   Usada por HookeJeeves.m; exemplo: ex03_2_5_hooke_jeeves.m
%   Implementação didática — algumas salvaguardas de software profissional
%   não estão incluídas.
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

if nargin < 5 || isempty(menor), menor = @(u, v) u < v; end   % a ordem habitual
xb = xb(:).';
n = numel(xb);
P = P(:).' .* ones(1, n);

xe = xb;  fe = fb;               % xe = cópia de xb
nfev = 0;
hist = zeros(0, n + 4);
ftrace = zeros(1, 0);

% --- núcleo (o dos slides, com as contagens; menor(u, v) = u < v) ---------
for j = 1:n                      % exploração (xe = cópia de xb)
  xp = xe;  xp(j) = xp(j) + P(j);  fp = f(xp);   % +P_j
  xm = xe;  xm(j) = xm(j) - P(j);  fm = f(xm);   % -P_j
  nfev = nfev + 2;  ftrace(end + 1:end + 2) = [fp(end), fm(end)];
  mais = menor(fp, fe) && ~menor(fm, fp);        % fp < fe e fp <= fm
  menos = ~mais && menor(fm, fe);                % fm < fe e fm < fp
  hist(end + 1, :) = [j, +1, xp, fp(end), mais];
  hist(end + 1, :) = [j, -1, xm, fm(end), menos];
  if mais                        % o melhor dos três: +P_j
    xe = xp;  fe = fp;
  elseif menos                   % o melhor dos três: -P_j
    xe = xm;  fe = fm;
  end
end
% ---------------------------------------------------------------------------

cols = cell(1, n + 4);
cols(1:2) = {'j', 'sinal'};
for j = 1:n, cols{2 + j} = sprintf('xp%d', j); end
cols(n + 3:n + 4) = {'f(xp)', 'aceite'};

info.nfev = nfev;
info.ngev = 0;
info.nhev = 0;
info.nit = n;
info.history = hist;
info.cols = cols;
info.ftrace = ftrace;
if menor(fe, fb)
  info.flag = 0;
  info.message = sprintf('a exploração melhorou: f = %.6g -> %.6g (%d avaliações)', fb(end), fe(end), nfev);
else
  info.flag = 1;
  info.message = sprintf('a exploração falhou: nenhuma das tentativas +-P_j melhorou (%d avaliações)', nfev);
end
end
