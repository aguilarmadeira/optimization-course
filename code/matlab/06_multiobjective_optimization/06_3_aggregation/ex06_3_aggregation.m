% EX06_3_AGGREGATION  Reproduz os exemplos do deck 6.3 (agregação de objetivos).
%
%   Exemplo convexo da UC (6.1):  f1 = x1^2 + x2^2, f2 = (x1 - 1)^2 + x2^2,
%   x em [0,1]^2, x1^2 + x2^2 <= 1.
%     1. «O que fazem os pesos?» e «Pesos uniformes ≠ pontos uniformes»:
%        soma ponderada com 11 pesos w1 = 0, 0.1, ..., 1, x0 = (0.5, 0.3),
%        arranque a quente; x1* = w2, f1 = w2^2, f2 = (1 - w2)^2 (tabela).
%   Exemplo não convexo da UC (6.2):  f1 = x1, f2 = 1 - x1^2 + x2,
%   x em [0,1] x [0, 0.6].
%     2. «A limitação fundamental»: com w1, w2 > 0 a soma ponderada só dá os
%        dois extremos (supported): 11, 100 e 10 000 pesos -> 2 pontos.
%     3. «Tchebycheff ponderado»: min theta s.a. w_i (f_i - z_i^id) <= theta,
%        z^id = (0, 0): encontra também os pontos non-supported.
%   4. «A escala importa»: projetos A, B, C; w = (0.9, 0.1) sem normalizar dá
%      81.4 / 85.8 / 100.5 (ganha A); normalizado com ideal e nadir, C / B / A.
%
%   Usa WeightedSum.m e WeightedTchebycheff.m (nesta pasta), com fmincon
%   (MATLAB) ou sqp (GNU Octave). O n_f = 75 do slide é do SLSQP do SciPy e só
%   se verifica na versão Python; aqui as contagens são as do solver usado e
%   imprimem-se só a título informativo (não entram na verificação).
%   Os slides usam vírgula decimal; aqui usa-se o ponto.
%
%   Otimização — deck 6.3.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

% v confere com o valor s do slide, mostrado com d casas decimais?
confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
perto = @(v, s, tol) all(abs(v(:) - s(:)) <= tol);
simnao = {'não', 'sim'};
limpa = @(v) v .* (abs(v) >= 5e-5) + 0;       % não imprimir -0.0000
ok = true;

fprintf('Otimização — deck 6.3: métodos de agregação de objetivos\n');

% ------------------------------------------------ 1. Exemplo convexo
F = @(x) [x(1)^2 + x(2)^2, (x(1) - 1)^2 + x(2)^2];
g = @(x) 1 - x(1)^2 - x(2)^2;                  % forma da UC: g >= 0
lb = [0 0];  ub = [1 1];  x0 = [0.5 0.3];
w1 = (0:0.1:1)';  W = [w1, 1 - w1];
fprintf('\n1. Exemplo convexo: min w1 f1 + w2 f2, 11 pesos, x0 = (%g, %g), arranque a quente\n', x0);
[x, fx, info] = WeightedSum(F, x0, lb, ub, g, W, true);
fprintf('   solver: %s (%s)\n', info.solver, info.message);
fprintf('   n_f = %d (%s; o slide dá n_f = 75 com o SLSQP do SciPy — só se verifica em Python)\n', ...
        info.nfev, info.solver);
T = [0.0 1.0 1.000 1.000 0.000;  0.1 0.9 0.900 0.810 0.010;  0.2 0.8 0.800 0.640 0.040
     0.3 0.7 0.700 0.490 0.090;  0.4 0.6 0.600 0.360 0.160;  0.5 0.5 0.500 0.250 0.250
     0.6 0.4 0.400 0.160 0.360;  0.7 0.3 0.300 0.090 0.490;  0.8 0.2 0.200 0.040 0.640
     0.9 0.1 0.100 0.010 0.810;  1.0 0.0 0.000 0.000 1.000];
c = [confere([W, x(:, 1), fx], T, 3), perto(x(:, 1), W(:, 2), 1e-5), perto(x(:, 2), 0, 1e-5), ...
     confere(x([11 6 1], :), [0 0; 0.5 0; 1 0], 3), info.flag == 0];
fprintf('  tabela (11 pesos): %s | x1* = w2: %s | x2* = 0: %s | w = (1;0), (0.5;0.5), (0;1) -> (0;0), (0.5;0), (1;0): %s | sem falhas: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);
d = diff(fx(:, 1));
fprintf('   pesos uniformes, pontos não uniformes: passos em f1 de %.2f (junto a f1 = 1) a %.2f (junto a f1 = 0)\n', ...
        -d(1), -d(end));

% ------------------------------------------------ 2. Exemplo não convexo
F = @(x) [x(1), 1 - x(1)^2 + x(2)];
lb = [0 0];  ub = [1 0.6];
X0 = [0.1 0.1; 0.5 0.1; 0.9 0.1];              % três pontos iniciais (solver local)
fprintf('\n2. Exemplo não convexo: soma ponderada, 3 pontos iniciais por peso, fica o melhor\n');
w1 = [0.2; 0.5; 0.8];
[x, fx] = WeightedSum(F, X0, lb, ub, [], [w1, 1 - w1]);
for k = 1:3
  fprintf('   w1 = %.1f: x* = (%.4f; %.4f), (f1, f2) = (%.4f; %.4f)\n', w1(k), limpa([x(k, :), fx(k, :)]));
end
[xa, fa] = WeightedSum(F, [0.5 0.1], lb, ub, [], [0.5 0.5]);
fprintf('   w1 = 0.5 só a partir de (0.5; 0.1): x = (%.4f; %.4f), f~ = %.4f', limpa(xa), [0.5 0.5]*fa');
if xa(1) > 0.01 && xa(1) < 0.99
  fprintf('\n   (cuidado: ponto estacionário de f~, côncava em x1 — é um máximo, não um mínimo)\n');
else
  fprintf(' (um extremo; o SLSQP do Python pára em (0.5; 0), máximo de f~ em x1)\n');
end
c = [];
for N = [11 100 10000]
  w1 = linspace(0, 1, N)';
  pos = w1 > 0 & w1 < 1;                       % w1, w2 > 0
  if N < 10000
    [x, fx] = WeightedSum(F, X0, lb, ub, [], [w1, 1 - w1]);
    como = 'solver';
  else                                         % 10 000 pesos: pesquisa exaustiva em x1
    t = linspace(0, 1, 1001);                  % (com w2 > 0 o ótimo tem x2 = 0)
    [~, j] = min(w1(pos)*t + (1 - w1(pos))*(1 - t.^2), [], 2);
    fx = zeros(N, 2);  fx(pos, :) = [t(j)', 1 - t(j)'.^2];
    como = 'grelha';
  end
  P = unique(round(fx(pos, :)*1e6)/1e6, 'rows');
  fprintf('   %5d pesos (%s), %5d com w1, w2 > 0 -> %d pontos distintos:%s\n', ...
          N, como, sum(pos), size(P, 1), sprintf(' (%g; %g)', P'));
  c(end + 1) = size(P, 1) == 2 && perto(sortrows(P), [0 1; 1 0], 1e-6);
end
fprintf('   (com w1 = 1, w2 = 0, x2 fica livre: o solver devolve (0; 1.1), só fracamente eficiente)\n');
fprintf('  11 pesos -> 2 pontos: %s | 100 -> 2: %s | 10 000 -> 2: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ 3. Tchebycheff ponderado
zid = [0 0];
w1 = (0:0.1:1)';  W = max([w1, 1 - w1], 1e-3);  % como no script das figuras
fprintf('\n3. Tchebycheff ponderado (min theta), z^id = (0, 0), 11 pesos (w_i >= 1e-3), a quente\n');
[x, fx, info] = WeightedTchebycheff(F, [0.5 0.3], lb, ub, [], W, zid, true);
fprintf('   solver: %s (%s); n_f = %d chamadas a F (informativo)\n', info.solver, info.message, info.nfev);
% exato: na frente f2 = 1 - f1^2, w1 f1 = w2 (1 - f1^2)
a = W(:, 1);  b = W(:, 2);
f1ex = (-a + sqrt(a.^2 + 4*b.^2))./(2*b);
nint = sum(fx(:, 1) > 0.01 & fx(:, 1) < 0.99);
fprintf('   f1 exato (w1 f1 = w2 (1 - f1^2)):%s\n', sprintf(' %.4f', f1ex));
fprintf('   %d dos 11 pontos são interiores (non-supported)\n', nint);
c = [perto(fx(:, 1), f1ex, 1e-3), perto(fx(:, 2), 1 - f1ex.^2, 2e-3), nint == 9, info.flag == 0];
fprintf('  pontos na frente (cruzamento w1 f1 = w2 f2): %s | f2 = 1 - f1^2: %s | 9 non-supported: %s | sem falhas: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ 4. A escala importa
fprintf('\n4. A escala importa: f1 = massa (kg), f2 = custo (€)\n');
Y = [1.5 800; 0.9 850; 0.5 1000];  nomes = 'ABC';
s = Y * [0.9; 0.1];
[~, j] = min(s);
fprintf('   w = (0.9; 0.1), sem normalizar: A = %.2f, B = %.2f, C = %.2f -> ganha %s\n', s, nomes(j));
zid = min(Y);  znad = max(Y);                  % os três são eficientes
Yn = (Y - zid) ./ (znad - zid);
Wn = [0.9 0.1; 0.5 0.5; 0.1 0.9];
[~, jn] = min(Yn * Wn', [], 1);
fprintf('   normalizado com z^id = (%g; %g), z^nad = (%g; %g):\n', zid, znad);
for k = 1:3
  v = Yn * Wn(k, :)';
  fprintf('     w = (%.1f; %.1f): A = %.3f, B = %.3f, C = %.3f\n', Wn(k, :), v);
end
fprintf('     escolhe %s / %s / %s\n', nomes(jn(1)), nomes(jn(2)), nomes(jn(3)));
c = [confere(s, [81.4; 85.8; 100.5], 1), j == 1, isequal(nomes(jn), 'CBA')];
fprintf('  81.4 / 85.8 / 100.5: %s | ganha A: %s | normalizado C / B / A: %s\n', simnao{c + 1});
ok = ok && all(c);

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
