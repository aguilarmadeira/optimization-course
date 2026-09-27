function hv = Hypervolume2D(F, r)
%HYPERVOLUME2D  Hipervolume em 2D (minimização): área dominada, limitada por r.
%
%   hv = Hypervolume2D(F, r)
%
%   Só contam os pontos que dominam estritamente r (f_1 < r_1 e f_2 < r_2);
%   retiram-se os dominados, ordena-se por f_1 e somam-se os retângulos
%       HV = soma_k (r_1 - f_1^(k)) * (f_2^(k-1) - f_2^(k)),  f_2^(0) = r_2.
%   Devolve 0 se nenhum ponto dominar r.
%
%   Entradas
%     F        matriz N x 2 dos objetivos
%     r        ponto de referência [r1 r2], pior do que todos
%
%   Otimização — deck 6.2 (frame «Medir a qualidade: o hipervolume») e
%   deck 6.6. Usada por NSGA2 e por ex06_6_nsga2.m.
%   Implementação didática — só para m = 2 objetivos.
%
%   J. F. A. Madeira — Licença MIT (ver LICENSE na raiz do repositório).

r = r(:)';
F = F(all(F < repmat(r, size(F, 1), 1), 2), :);   % só os que dominam r
hv = 0;
if isempty(F), return, end
N = size(F, 1);
nd = true(N, 1);
for p = 1:N                                        % retirar os dominados
  Fp = repmat(F(p, :), N, 1);
  nd(p) = ~any(all(F <= Fp, 2) & any(F < Fp, 2));
end
F = F(nd, :);
[~, o] = sort(F(:, 1));
F = F(o, :);
prev = r(2);
for k = 1:size(F, 1)
  hv = hv + (r(1) - F(k, 1))*(prev - F(k, 2));
  prev = F(k, 2);
end
end
