% EX03_3_3_CONJUGATE_GRADIENTS  Reproduz os exemplos do deck 3.3.3
% (gradientes conjugados, Fletcher-Reeves).
%
%   Quadrática f(x) = (x1^2 + 9 x2^2)/2 de x0 = (9, 1), pesquisa em linha
%       exata: o 1.º passo é o do gradiente, x1 = (7.2; -0.8);
%       beta_1 = 103.68/162 = 0.64, d1 = (-12.96; 1.44), d0'*A*d1 = 0
%       (A-conjugadas), alpha_1 = 5/9, x2 = (0, 0): 2 iterações, n_g = 3.
%   Rosenbrock de (-1.5; 2), ||g|| < 1e-4: com reinício a cada n = 2, 36 it.
%       (n_g = 37); sem reinício, 102 it.
%   «No computador»: com uma pesquisa do tipo fminbnd de tolerância por
%       omissão, FR no Rosenbrock pode não convergir (informativo).
%
%   Imprime as tabelas, os valores finais e as contagens, e no fim compara
%   com os slides (vírgula decimal nos slides, ponto aqui).
%
%   Otimização — deck 3.3.3.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
simnao = {'não', 'sim'};
ok = true;

fprintf('Otimização — deck 3.3.3: gradientes conjugados (Fletcher-Reeves)\n');

% ------------------------------------------------ quadrática
A = diag([1 9]);
quad = @(x) 0.5*(x(1)^2 + 9*x(2)^2);          % = x'*A*x/2
gq = @(x) [x(1); 9*x(2)];
x0 = [9; 1];
fprintf('\nQuadrática (x1^2 + 9 x2^2)/2 de (9, 1), tolg = 1e-6, pesquisa em linha de Brent (tol 1e-10):\n');
[x, ~, info] = ConjGrad(quad, gq, x0, 1e-6, 100, struct('verbose', true));
H = info.history;                              % [k x1 x2 f ||g|| alpha beta d1 d2]
g0 = gq(H(1, 2:3)');  g1 = gq(H(2, 2:3)');
d0 = H(1, 8:9)';  d1 = H(2, 8:9)';
fprintf('||g1||^2/||g0||^2 = %.2f/%.0f = %.4f; d0''*A*d1 = %.1e; alpha_1 = %.6f (5/9 = %.6f)\n', ...
        g1'*g1, g0'*g0, H(2, 7), d0'*A*d1, H(3, 6), 5/9);
fprintf('%s; n_g = %d, n_f = %d\n', info.message, info.ngev, info.nfev);
c = [confere(H(2, 2:3), [7.2 -0.8], 4), confere(g1'*g1, 103.68, 2), confere(H(2, 7), 0.64, 4), ...
     confere(d1, [-12.96; 1.44], 2), abs(d0'*A*d1) < 1e-4, confere(H(3, 6), 5/9, 6), ...
     confere(x, [0; 0], 6), info.nit == 2, info.ngev == 3];
fprintf(['  x1: %s | ||g1||^2 = 103.68: %s | beta_1 = 0.64: %s | d1: %s | d0''Ad1 = 0: %s | alpha_1 = 5/9: %s', ...
         ' | x2 = (0, 0): %s | 2 it.: %s | n_g = 3: %s\n'], simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ Rosenbrock
rosen = @(x) (1 - x(1))^2 + 100*(x(2) - x(1)^2)^2;
grosen = @(x) [-2*(1 - x(1)) - 400*x(1)*(x(2) - x(1)^2); 200*(x(2) - x(1)^2)];
xr = [-1.5; 2];
fprintf('\nRosenbrock de (-1.5; 2), tolg = 1e-4:\n');
[~, fc, ic] = ConjGrad(rosen, grosen, xr, 1e-4, 500);                        % reinício a cada n = 2
[~, fs, is] = ConjGrad(rosen, grosen, xr, 1e-4, 500, struct('restart', 0));  % sem reinício periódico
fprintf('FR com reinício a cada n = 2: %d it., n_g = %d, n_f = %d, f = %.1e (%d reinícios)\n', ...
        ic.nit, ic.ngev, ic.nfev, fc, ic.nrestart);
fprintf('FR sem reinício periódico:    %d it., n_g = %d, n_f = %d, f = %.1e (%d reinícios por g''*d >= 0)\n', ...
        is.nit, is.ngev, is.nfev, fs, is.nrestart);
c = [ic.nit == 36, ic.ngev == 37, is.nit == 102, ic.flag == 0 && is.flag == 0];
fprintf('  36 it. com reinício: %s | n_g = 37: %s | 102 it. sem: %s | convergem: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------- pesquisa tipo fminbnd (informativo)
[~, fb, ib] = ConjGrad(rosen, grosen, xr, 1e-4, 500, struct('ls', 'fminbnd'));
fprintf('\nCom a pesquisa do tipo fminbnd (TolX = 1e-4) e reinício a cada n: %s; f = %.2g\n', ib.message, fb);
if ib.flag, s = 'não convergiu'; else, s = 'convergiu'; end
fprintf('  (o slide diz que «pode não convergir»: %s)\n', s);

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
