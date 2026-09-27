function [xe, fe, info] = Exploratory(f, xb, fb, P)
%EXPLORATORY  Movimento exploratório de Hooke–Jeeves em torno de xb.
%
%   [xe, fe, info] = Exploratory(f, xb, fb, P)
%
%   Percorre as coordenadas j = 1, ..., n, sequencialmente: tenta x_j + P_j;
%   se não melhora, tenta x_j - P_j; se nenhuma melhora, deixa x_j como estava.
%   Cada teste parte do melhor ponto encontrado até aí.
%
%   Entradas
%     f     função (handle) de um vetor linha x (1 x n)
%     xb    ponto base (vetor; é tratado como linha)
%     fb    f(xb), já conhecido (NÃO é avaliado aqui nem contado)
%     P     perturbações, uma por coordenada (escalar = igual em todas)
%
%   Saídas
%     xe    ponto explorado (vetor linha); xe = xb se a exploração falhou
%     fe    f(xe)
%     info  estrutura com
%       .nfev     avaliações de f: entre n (todas as tentativas +P_j funcionam)
%                 e 2n (todas falham)
%       .ngev, .nhev   0
%       .nit      coordenadas percorridas (= n)
%       .history  uma linha por avaliação: [j  sinal  xp(1:n)  f(xp)  aceite]
%                 sinal = +1 ou -1 (x_j + P_j ou x_j - P_j), aceite = 1 ou 0
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

xb = xb(:).';
n = numel(xb);
P = P(:).' .* ones(1, n);

xe = xb;  fe = fb;               % xe = cópia de xb
nfev = 0;
hist = zeros(0, n + 4);
ftrace = zeros(1, 0);

% --- núcleo (igual ao dos slides, com as contagens) ------------------------
for j = 1:n                      % exploração (xe = cópia de xb)
  xp = xe;  xp(j) = xp(j) + P(j);  fp = f(xp);   % xp: teste
  nfev = nfev + 1;  hist(end + 1, :) = [j, +1, xp, fp, fp < fe];  ftrace(end + 1) = fp;
  if fp >= fe
    xp(j) = xe(j) - P(j);  fp = f(xp);
    nfev = nfev + 1;  hist(end + 1, :) = [j, -1, xp, fp, fp < fe];  ftrace(end + 1) = fp;
  end
  if fp < fe, xe = xp;  fe = fp;  end
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
if fe < fb
  info.flag = 0;
  info.message = sprintf('a exploração melhorou: f = %.6g -> %.6g (%d avaliações)', fb, fe, nfev);
else
  info.flag = 1;
  info.message = sprintf('a exploração falhou: nenhuma das tentativas +-P_j melhorou (%d avaliações)', nfev);
end
end
