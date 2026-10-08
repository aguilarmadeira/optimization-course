% EX03_2_3_NELDER_MEAD  Reproduz os exemplos do deck 3.2.3 (simplex de Nelder–Mead).
%
%   Notação e paragem das aulas: xl, xg, xh; gamma = 2, beta = 0.5;
%   análise do erro Q = sqrt(sum (f(x_i) - f(xc))^2 / (n+1)), parar se Q <= eps.
%   1. Uma iteração à mão: f = x1^2 + x2^2, simplex (-1,0), (1,2), (-1,3):
%      reflexão aceite, Q = 2.380, 2 avaliações (xc e xr).
%   2. Himmelblau a partir de (0,0), simplex (0,0), (1.2,0), (0,1.2), eps = 1e-6:
%      tabela k = 0..12 (com Q); 31 iterações, n_f = 92, x* = (3,2).
%      Com o simplex (0,0), (-1.2,0), (0,-1.2) chega-se a (-2.81; 3.13).
%   3. Rosenbrock a partir de (-1.5, 2): n_f até f < 1e-4 = 244 (simplex de
%      arestas 0.5) e 211 (arestas de 5 %); sem as avaliações de f(xc), isto é,
%      com a paragem do Nelder–Mead padrão (tolx, tolf, como o fminsearch):
%      166 e 139.
%
%   O deck mostra também fminsearch; aqui usa-se a implementação da UC,
%   NelderMead, com o algoritmo do deck (fminsearch tem outro simplex inicial
%   e outro critério de paragem: mesma ideia, contagens diferentes).
%   Os slides usam vírgula decimal; aqui usa-se o ponto.
%
%   Otimização — deck 3.2.3.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

% v confere com o valor s do slide, mostrado com d casas decimais?
confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
simnao = {'não', 'sim'};
ok = true;
aulas = struct('eps', 1e-6);          % gamma = 2 e beta = 0.5 por omissão

fprintf('Otimização — deck 3.2.3: método do simplex de Nelder–Mead\n');

sq     = @(x) x(1)^2 + x(2)^2;
himmel = @(x) (x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2;
ros    = @(x) (1 - x(1))^2 + 100*(x(2) - x(1)^2)^2;

% ------------------------------------------------------------ 1. à mão
fprintf('\n1. Uma iteração à mão: f = x1^2 + x2^2, simplex (-1,0), (1,2), (-1,3)\n');
[~, ~, info] = NelderMead(sq, [-1 0; 1 2; -1 3], [], [], 1, [], true, struct('eps', 0));
S = info.simplex;  Fs = [sq(S(1, :)), sq(S(2, :)), sq(S(3, :))];
fprintf('novo simplex: (%g,%g), (%g,%g), (%g,%g) com f = %g, %g, %g;  Q = %.3f;  avaliações nesta iteração: %d\n', ...
        S', Fs, info.Q, info.nfev - 3);
c = [strcmp(info.ops{2}, 'reflexão'), confere(S, [-1 0; 1 -1; 1 2], 4), ...
     confere(Fs, [1 2 5], 4), confere(info.Q, 2.380, 3), info.nfev - 3 == 2];
fprintf('  reflexão aceite: %s | novo simplex: %s | f = 1, 2, 5: %s | Q = 2.380: %s | 2 avaliações: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------------ 2. Himmelblau
fprintf('\n2. Himmelblau a partir de (0,0), simplex (0,0), (1.2,0), (0,1.2), eps = 1e-6\n');
[x, fx, info] = NelderMead(himmel, [0 0; 1.2 0; 0 1.2], [], [], [], [], false, aulas);
fprintf('primeiras 13 linhas de info.history [k, n_f, f(x_l), f(x_h), Q, x_l] e operação:\n');
fprintf('%4s %5s %12s %12s %10s %10s %10s  operação\n', info.cols{:});
for i = 1:13
  fprintf('%4d %5d %12.4f %12.4f %10.4f %10.4f %10.4f  %s\n', info.history(i, :), info.ops{i});
end
fprintf('... %d iterações, n_f = %d, x* = (%.6f, %.6f), f = %.1e  (%s)\n', ...
        info.nit, info.nfev, x(1), x(2), fx, info.message);
% tabela do slide (k = 0..12): operação, x_l, f(x_l), f(x_h), Q
Tops = {'inicial', 'expansão', 'reflexão', 'reflexão', 'contr. int.', 'contr. int.', 'contr. ext.', ...
        'contr. int.', 'contr. int.', 'contr. int.', 'contr. int.', 'contr. ext.', 'contr. int.'};
Txl = [1.200 0.000; 1.800 1.800; 3.000 0.600; 3.000 0.600; 2.550 1.650; 3.188 1.762; 3.188 1.762
       3.188 1.762; 2.884 1.921; 2.919 2.050; 2.919 2.050; 3.031 1.983; 3.031 1.983];
Tfl = [125.0336 39.3632 15.2096 15.2096 11.0925 1.3499 1.3499 1.3499 0.7628 0.1972 0.1972 0.0303 0.0303];
Tfh = [170.000 126.954 125.034 39.363 24.579 15.210 11.093 2.965 1.604 1.350 0.763 0.217 0.197];
TQ = [57.286 52.768 9.943 14.104 7.211 4.981 2.095 1.286 0.643 0.262 0.130 0.104];
H = info.history(1:13, :);
c = [isequal(info.ops(1:13), Tops), confere(H(:, 6:7), Txl, 3), confere(H(:, 3), Tfl, 4), ...
     confere(H(:, 4), Tfh, 3), confere(H(2:13, 5), TQ, 3), info.nit == 31, info.nfev == 92, ...
     confere(x, [3 2], 3), fx < 1e-5];
fprintf(['  tabela k = 0..12: operação: %s | x_l: %s | f(x_l): %s | f(x_h): %s | Q: %s\n' ...
         '  31 it.: %s | n_f = 92: %s | x* = (3,2): %s | f ~ 0: %s\n'], simnao{c + 1});
ok = ok && all(c);

[x2, f2] = NelderMead(himmel, [0 0; -1.2 0; 0 -1.2], [], [], [], [], false, aulas);
fprintf('simplex (0,0), (-1.2,0), (0,-1.2): x* = (%.2f; %.2f), f = %.1e\n', x2(1), x2(2), f2);
c = confere(x2, [-2.81 3.13], 2);
fprintf('  chega a (-2.81; 3.13): %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------------------ 3. Rosenbrock
x0 = [-1.5 2];
fprintf('\n3. Rosenbrock a partir de (-1.5, 2): n_f até f < 1e-4 (info.nhit)\n');
[xa, ~, ia] = NelderMead(ros, [x0; x0 + [0.5 0]; x0 + [0 0.5]], [], [], [], 1e-4, false, aulas);
[xb, ~, ib] = NelderMead(ros, x0, [], [], [], 1e-4, false, aulas);   % simplex por omissão: h_i = 0.05 max(1,|x0_i|)
fprintf('simplex de arestas 0.5:     n_f até f < 1e-4 = %d  (final: %d it., n_f = %d, x = (%.4f, %.4f))\n', ...
        ia.nhit, ia.nit, ia.nfev, xa);
fprintf('simplex de arestas de 5 %%:  n_f até f < 1e-4 = %d  (final: %d it., n_f = %d, x = (%.4f, %.4f))\n', ...
        ib.nhit, ib.nit, ib.nfev, xb);
[~, ~, pa] = NelderMead(ros, [x0; x0 + [0.5 0]; x0 + [0 0.5]], [], [], [], 1e-4);   % paragem padrão, sem f(xc)
[~, ~, pb] = NelderMead(ros, x0, [], [], [], 1e-4);
fprintf('sem as avaliações de f(x_c) (paragem do Nelder–Mead padrão): %d e %d\n', pa.nhit, pb.nhit);
c = [ia.nhit == 244, ib.nhit == 211, pa.nhit == 166, pb.nhit == 139];
fprintf('  n_f = 244: %s | n_f = 211: %s | sem f(x_c): 166: %s, 139: %s\n', simnao{c + 1});
ok = ok && all(c);

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
