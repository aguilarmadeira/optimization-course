% EX03_3_1_STEEPEST_DESCENT  Reproduz os exemplos do deck 3.3.1 (método do
% gradiente) e o do gradiente numérico do deck 3.3.
%
%   «Uma iteração à mão»: f(x) = (x1^2 + 9 x2^2)/2, x0 = (9, 1): alpha_0 = 0.2,
%       x1 = (7.2; -0.8), f = 28.8, g1 = (7.2; -7.2), g1'*d0 = 0.
%   «O percurso em ziguezague»: tabela k = 0..8; 74 iterações para
%       ||g|| < 1e-6, n_g = 75; f multiplica por 0.64 em cada iteração.
%   «No computador»: pesquisa em linha de alta precisão (a das figuras, cerca
%       de 24 avaliações de f por pesquisa) e do tipo fminbnd (6 por
%       pesquisa): as mesmas 74 iterações.
%   Deck 3.3: gradiente de Rosenbrock em (-1.5; 2) = (-155, -50); com
%       h = 1e-4 o erro das diferenças centrais é 6e-6.
%   «Rosenbrock»: de (-1.5; 2), 3000 iterações, ainda f = 2.4e-4,
%       x ~ (0.98; 0.97).
%
%   Imprime as tabelas, os valores finais e as contagens, e no fim compara
%   com os slides (vírgula decimal nos slides, ponto aqui).
%
%   Otimização — deck 3.3.1.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
simnao = {'não', 'sim'};
ok = true;

fprintf('Otimização — deck 3.3.1: método do gradiente\n');

f = @(x) 0.5*(x(1)^2 + 9*x(2)^2);        % = x'*A*x/2, A = diag(1, 9)
grad = @(x) [x(1); 9*x(2)];
x0 = [9; 1];

% ------------------------------------------------ uma iteração à mão
g0 = grad(x0);  d0 = -g0;
[x1, ~, i1] = SteepestDescent(f, grad, x0, 1e-6, 1);    % uma só iteração
g1 = grad(x1);  alpha0 = i1.history(2, 6);
fprintf('\nUma iteração à mão: g0 = (%g, %g), d0 = (%g, %g)\n', g0, d0);
fprintf('alpha_0 = %.4f (exato: 162/810 = 0.2), x1 = (%.4f; %.4f), f(x1) = %.4f (era %g)\n', ...
        alpha0, x1, f(x1), f(x0));
fprintf('g1 = (%.4f; %.4f), g1''*d0 = %.1e (ortogonais)\n', g1, g1'*d0);
c = [confere(alpha0, 0.2, 4), confere(x1, [7.2; -0.8], 4), confere(f(x1), 28.8, 4), ...
     confere(g1, [7.2; -7.2], 4), abs(g1'*d0) < 1e-6];
fprintf('  alpha_0 = 0.2: %s | x1: %s | f = 28.8: %s | g1: %s | g1''*d0 = 0: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ percurso em ziguezague
tolg = 1e-6;
fprintf('\nZiguezague: f(x) = (x1^2 + 9 x2^2)/2, x0 = (9, 1), tolg = %g, pesquisa em linha de Brent (tol 1e-10)\n', tolg);
[x, fx, info] = SteepestDescent(f, grad, x0, tolg, 1000);
H = info.history;
fprintf('%4s %22s %12s %10s %9s\n', 'k', 'x_k', 'f(x_k)', '||g_k||', 'alpha_k-1');
for r = 1:9
  fprintf('%4d  (%9.4f; %9.4f) %12.5f %10.4f %9.4f\n', H(r, :));
end
fprintf('...\n');
fprintf('%4d  (%9.2e; %9.2e) %12.2e %10.2e %9.4f\n', H(end, :));
razao = H(2:end, 4)./H(1:end-1, 4);
fprintf('%s: %d iterações, n_g = %d, n_f = %d (%d nas pesquisas em linha + %d iteradas)\n', ...
        info.message, info.nit, info.ngev, info.nfev, info.nfev_ls, info.nit + 1);
fprintf('média de %.1f avaliações de f por pesquisa em linha (o slide diz «cerca de 24»)\n', info.nfev_ls/info.nit);
fprintf('f_{k+1}/f_k: mín. %.4f, máx. %.4f (cota de Kantorovich ((9-1)/(9+1))^2 = 0.64)\n', min(razao), max(razao));
T = [9.0000  1.0000 45.00000 12.7279 NaN
     7.2000 -0.8000 28.80000 10.1823 0.2000
     5.7600  0.6400 18.43200  8.1459 0.2000
     4.6080 -0.5120 11.79648  6.5167 0.2000
     3.6864  0.4096  7.54975  5.2134 0.2000
     2.9491 -0.3277  4.83184  4.1707 0.2000
     2.3593  0.2621  3.09238  3.3365 0.2000
     1.8874 -0.2097  1.97912  2.6692 0.2000
     1.5099  0.1678  1.26664  2.1354 0.2000];
c = [confere(H(1:9, 2:3), T(:, 1:2), 4) && confere(H(1:9, 4), T(:, 3), 5) && confere(H(1:9, 5), T(:, 4), 4) ...
       && confere(H(2:9, 6), T(2:end, 5), 4), ...
     info.nit == 74, info.ngev == 75, info.flag == 0, confere([min(razao) max(razao)], [0.64 0.64], 4), ...
     round(info.nfev_ls/info.nit) == 24];
fprintf('  tabela k = 0..8: %s | 74 it.: %s | n_g = 75: %s | ||g|| < 1e-6: %s | razão 0.64: %s | cerca de 24 por pesquisa: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------ pesquisa em linha tipo fminbnd
opts.ls = 'fminbnd';
[~, ~, ib] = SteepestDescent(f, grad, x0, tolg, 1000, opts);
fprintf('\nCom a pesquisa do tipo fminbnd em [0, 1] (TolX = 1e-4): %d iterações, n_g = %d, n_f = %d; %.1f avaliações por pesquisa\n', ...
        ib.nit, ib.ngev, ib.nfev, ib.nfev_ls/ib.nit);
c = [ib.nit == 74, ib.ngev == 75, ib.nfev_ls/ib.nit == 6];
fprintf('  74 it.: %s | n_g = 75: %s | 6 avaliações por pesquisa: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ gradiente numérico
[~, ~, in] = SteepestDescent(f, [], x0, tolg, 1000);
fprintf('\nCom gradiente numérico (grad = []: GradFD, 2n = 4 avaliações de f por gradiente):\n');
fprintf('%d iterações, n_g = %d, n_f = %d = %d (pesquisas) + %d (iteradas) + 4 x %d (gradientes)\n', ...
        in.nit, in.ngev, in.nfev, in.nfev_ls, in.nit + 1, in.nit + 1);
c = [in.nit == 74, in.ngev == 0, in.nfev == in.nfev_ls + (in.nit + 1) + 4*(in.nit + 1)];
fprintf('  74 it.: %s | n_g = 0: %s | 2n por gradiente em n_f: %s\n', simnao{c + 1});
ok = ok && all(c);

rosen = @(x) (1 - x(1))^2 + 100*(x(2) - x(1)^2)^2;
grosen = @(x) [-2*(1 - x(1)) - 400*x(1)*(x(2) - x(1)^2); 200*(x(2) - x(1)^2)];
xr = [-1.5; 2];
ge = grosen(xr);
fprintf('\nDeck 3.3 — gradiente de Rosenbrock em (-1.5; 2): exato (%g, %g)\n', ge);
hs = [1e-2 1e-4 1e-6];  erros = zeros(1, 3);
for j = 1:3
  [gn, nf] = GradFD(rosen, xr, hs(j));
  erros(j) = norm(gn - ge);
  fprintf('  h = %g: (%.6f, %.6f), erro %.1e, %d avaliações de f\n', hs(j), gn, erros(j), nf);
end
gd = GradFD(rosen, xr);
fprintf('  h por omissão (max(0.01|x_i|, 1e-4)): erro %.1e\n', norm(gd - ge));
c = [confere(ge, [-155; -50], 0), confere(erros(2), 6e-6, 6), nf == 4];
fprintf('  grad = (-155, -50): %s | erro 6e-6 com h = 1e-4: %s | 2n = 4 avaliações: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ Rosenbrock
fprintf('\nRosenbrock de (-1.5; 2), tolg = 1e-4, kmax = 3000:\n');
[xR, fR, iR] = SteepestDescent(rosen, grosen, xr, 1e-4, 3000);
fprintf('%s\nx = (%.4f; %.4f), f = %.2e, ||g|| no último gradiente = %.2e; n_g = %d, n_f = %d\n', ...
        iR.message, xR, fR, iR.history(end - 1, 5), iR.ngev, iR.nfev);
c = [iR.nit == 3000, iR.flag == 1, confere(fR, 2.4e-4, 5), confere(xR, [0.98; 0.97], 2)];
fprintf('  3000 it.: %s | não converge: %s | f = 2.4e-4: %s | x ~ (0.98; 0.97): %s\n', simnao{c + 1});
ok = ok && all(c);

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
