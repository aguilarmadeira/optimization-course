% EX03_3_2_NEWTON  Reproduz os exemplos do deck 3.3.2 (método de Newton).
%
%   «Podemos usar a curvatura?»: f(x) = (x1^2 + 9 x2^2)/2, x0 = (9, 1):
%       H d = -g dá d0 = (-9, -1) e x1 = (0, 0) -- uma iteração.
%   «Rosenbrock: muito rápido, mas não monótono»: Newton puro de (-1.5; 2),
%       tabela k = 0..7 (f sobe para 751 e para 40), 7 iterações, f < 1e-16.
%   «Newton procura g = 0»: Himmelblau de (0, 0), lambda(H0) = (-42; -26);
%       o puro vai em 4 iterações para o máximo local (-0.27; -0.92),
%       f = 181.6; o amortecido deteta g'*d >= 0, usa -g e chega a (3, 2)
%       em 5 iterações.
%   «Newton amortecido»: Rosenbrock, 14 iterações (monótono); os últimos
%       minimizantes da linha são 0.91; 0.98; 1.00.
%   Deck 3.3.3 (quadro de consulta): amortecido com ||g|| < 1e-4: 14 it.,
%       n_g = 15, n_H = 14.
%
%   Paragem ||g|| <= tolg = 1e-8 (a das figuras), salvo indicação.
%   Imprime as tabelas, os valores finais e as contagens, e no fim compara
%   com os slides (vírgula decimal nos slides, ponto aqui).
%
%   Otimização — deck 3.3.2.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
confere_sig = @(v, s) all(abs(v(:) - s(:)) <= 0.5*10.^(floor(log10(abs(s(:)))) - 2)*(1 + 1e-6));  % 3 alg. sig.
simnao = {'não', 'sim'};
nomes = {'Newton', 'Newton com H + mu I', '-g (salvaguarda)'};
ok = true;

fprintf('Otimização — deck 3.3.2: método de Newton\n');

% ------------------------------------------------ quadrática
A = diag([1 9]);
quad = @(x) 0.5*(x(1)^2 + 9*x(2)^2);          % = x'*A*x/2
gq = @(x) [x(1); 9*x(2)];
Hq = @(x) A;
x0 = [9; 1];
d0 = -(A \ gq(x0));
[x, ~, info] = NewtonND(quad, gq, Hq, x0, 1e-8, 50);
fprintf('\nQuadrática (x1^2 + 9 x2^2)/2 de (9, 1): gradiente -g0 = (%g, %g); Newton H d = -g: d0 = (%g, %g)\n', ...
        -gq(x0), d0);
fprintf('x1 = (%g, %g): %d iteração; n_g = %d, n_H = %d\n', x, info.nit, info.ngev, info.nhev);
c = [all(abs(d0 - [-9; -1]) < 1e-12), all(abs(x) < 1e-12), info.nit == 1];
fprintf('  d0 = (-9, -1): %s | x1 = (0, 0): %s | 1 it.: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ Rosenbrock, puro
rosen = @(x) (1 - x(1))^2 + 100*(x(2) - x(1)^2)^2;
grosen = @(x) [-2*(1 - x(1)) - 400*x(1)*(x(2) - x(1)^2); 200*(x(2) - x(1)^2)];
Hrosen = @(x) [2 - 400*(x(2) - 3*x(1)^2), -400*x(1); -400*x(1), 200];
xr = [-1.5; 2];
fprintf('\nRosenbrock de (-1.5; 2), Newton puro (tolg = 1e-8):\n');
opts = struct('verbose', true);
[x, fx, info] = NewtonND(rosen, grosen, Hrosen, xr, 1e-8, 50, opts);
H = info.history;
fprintf('%s; n_g = %d, n_H = %d, n_f = %d (f só nos iterandos)\n', info.message, info.ngev, info.nhev, info.nfev);
ks = find(diff(H(:, 4)) > 0);                 % f sobe de k-1 para k
sub = arrayfun(@(j) sprintf('k = %d (f = %.2f)', j, H(j + 1, 4)), ks(:)', 'UniformOutput', false);
fprintf('f sobe em: %s\n', strjoin(sub, ', '));
Tx = [-1.5000 2.0000; -1.4510 2.1029; 0.2044 -2.6986; 0.2059 0.0424
       0.9997 0.3692;  0.9997 0.9993; 1.0000 1.0000; 1.0000 1.0000];
Tf = [12.500000; 6.007882; 751.610005; 0.630622; 39.701728];
Tf56 = [1.09e-7; 1.20e-12];
Tg = [1.63e2; 6.31e0; 5.92e2; 1.59e0; 2.82e2; 6.61e-4; 4.89e-5];
c = [size(H, 1) == 8 && confere(H(:, 2:3), Tx, 4) && confere(H(1:5, 4), Tf, 6) ...
       && confere_sig(H(6:7, 4), Tf56) && H(8, 4) == 0 && confere_sig(H(1:7, 5), Tg) && H(8, 5) == 0, ...
     info.nit == 7, fx < 1e-16, isequal(ks(:)', [2 4])];
fprintf('  tabela k = 0..7: %s | 7 it.: %s | f < 1e-16: %s | não monótono (sobe em k = 2 e 4): %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ Rosenbrock, amortecido
fprintf('\nRosenbrock de (-1.5; 2), Newton amortecido (salvaguarda g''*d >= 0 -> -g; pesquisa em linha de Brent):\n');
opts = struct('damped', true, 'verbose', true);
[x, fx, info] = NewtonND(rosen, grosen, Hrosen, xr, 1e-8, 50, opts);
H = info.history;
mono = all(diff(H(:, 4)) <= 0);
fprintf('%s; n_g = %d, n_H = %d, n_f = %d (%d nas pesquisas em linha)\n', ...
        info.message, info.ngev, info.nhev, info.nfev, info.nfev_ls);
fprintf('monótono: %s; últimos minimizantes da linha: %.2f; %.2f; %.2f\n', simnao{mono + 1}, H(end-2:end, 6));
c = [info.nit == 14, mono, confere(H(end-2:end, 6), [0.91; 0.98; 1.00], 2)];
fprintf('  14 it.: %s | monótono: %s | alpha = 0.91; 0.98; 1.00: %s\n', simnao{c + 1});
ok = ok && all(c);
opts = struct('damped', true);
[~, ~, i4] = NewtonND(rosen, grosen, Hrosen, xr, 1e-4, 50, opts);
fprintf('Com tolg = 1e-4 (quadro de consulta do 3.3.3): %d it., n_g = %d, n_H = %d\n', i4.nit, i4.ngev, i4.nhev);
c = [i4.nit == 14, i4.ngev == 15, i4.nhev == 14];
fprintf('  14 it.: %s | n_g = 15: %s | n_H = 14: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ Himmelblau
himmel = @(x) (x(1)^2 + x(2) - 11)^2 + (x(1) + x(2)^2 - 7)^2;
ghimmel = @(x) [4*x(1)*(x(1)^2 + x(2) - 11) + 2*(x(1) + x(2)^2 - 7);
                2*(x(1)^2 + x(2) - 11) + 4*x(2)*(x(1) + x(2)^2 - 7)];
Hhimmel = @(x) [12*x(1)^2 + 4*x(2) - 42, 4*x(1) + 4*x(2);
                4*x(1) + 4*x(2), 4*x(1) + 12*x(2)^2 - 26];
xh = [0; 0];
lam0 = eig(Hhimmel(xh));
fprintf('\nHimmelblau de (0, 0): lambda(H0) = (%g; %g) -- definida negativa\n', lam0);
[xp, fp, ip] = NewtonND(himmel, ghimmel, Hhimmel, xh, 1e-8, 50);
lamx = eig(Hhimmel(xp));
fprintf('Puro: %d it. -> x = (%.4f; %.4f), f = %.4f, lambda(H) = (%.1f; %.1f) < 0: máximo local\n', ...
        ip.nit, xp, fp, lamx);
[xd, fd, id] = NewtonND(himmel, ghimmel, Hhimmel, xh, 1e-8, 50, struct('damped', true));
fprintf('Amortecido: direção na 1.ª iteração: %s; %d it. -> x = (%.4f; %.4f), f = %.1e\n', ...
        nomes{id.history(2, 8) + 1}, id.nit, xd, fd);
c = [all(abs(sort(lam0) - [-42; -26]) < 1e-12), ip.nit == 4, confere(xp, [-0.27; -0.92], 2), ...
     confere(fp, 181.6, 1), all(lamx < 0), id.history(2, 8) == 2, id.nit == 5, all(abs(xd - [3; 2]) < 1e-8)];
fprintf(['  lambda(H0) = (-42; -26): %s | puro 4 it.: %s | (-0.27; -0.92): %s | f = 181.6: %s | máximo: %s', ...
         ' | amortecido usa -g: %s | 5 it.: %s | (3, 2): %s\n'], simnao{c + 1});
ok = ok && all(c);

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
