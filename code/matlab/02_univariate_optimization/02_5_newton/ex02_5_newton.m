% EX02_5_NEWTON  Reproduz os exemplos do deck 2.5 (Newton-Raphson).
%
%   Exemplo-guia: f(x) = x^2 + 54/x, x0 = 2.5: x1 = 2.908; 4 iterações até
%       x* = 3 (erro 2.4e-12); n_g = n_H = 1 por iteração (+ f'(x0)).
%       Paragem |f'| < tolg com tolg = 1e-8 (o slide não fixa tolg; qualquer
%       valor entre 1.5e-11 e 1.6e-5 dá as mesmas 4 iterações).
%   «Convergência quadrática»: x^3 - 2x^2 + 4, x0 = 1, tolg = 1e-3 (3 it.).
%   «Newton procura f' = 0»: maximizar 8 + 5x - 3x^4 - 2x^6, x0 = 0.5,
%       tolg = 1e-7 (5 it.).
%   «Os quatro métodos»: x^2 + 54/x a partir de x0 = 1 (erro após 5 e 7 it.).
%
%   Imprime as tabelas como nos slides, os valores finais e as contagens, e
%   no fim compara com os slides (vírgula decimal nos slides, ponto aqui).
%
%   Otimização — deck 2.5.  J. F. A. Madeira — Licença MIT.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..')); uc_setup;   % caminhos do código da UC

confere = @(v, s, d) all(abs(v(:) - s(:)) <= 0.5*10.^(-d(:))*(1 + 1e-6));
simnao = {'não', 'sim'};
ok = true;

fprintf('Otimização — deck 2.5: Newton-Raphson\n');

% ------------------------------------------------------- exemplo-guia
f   = @(x) x.^2 + 54./x;
df  = @(x) 2*x - 54./x.^2;
ddf = @(x) 2 + 108./x.^3;
x0 = 2.5;  tolg = 1e-8;  kmax = 20;
fprintf('\nExemplo-guia: f(x) = x^2 + 54/x, x0 = %g, tolg = %g\n', x0, tolg);
[x, fx, info] = Newton1D(f, df, ddf, x0, tolg, kmax);
H = info.history;
fprintf('%3s %11s %11s %9s %11s %12s\n', 'k', 'x_k', 'df(x_k)', 'ddf(x_k)', 'x_k+1', '|x_k+1 - 3|');
for r = 1:size(H, 1)
  fprintf('%3d %11.7f %+11.5f %9.4f %11.7f %12.1e\n', H(r, [1 2 4 5 6]), abs(H(r, 6) - 3));
end
fprintf('x* = %.7f, f(x*) = %.4f; %d iterações, n_g = %d, n_H = %d, n_f = %d\n', ...
        x, fx, info.nit, info.ngev, info.nhev, info.nfev);
[x5, ~, info5] = Newton1D([], df, ddf, x0, 1e-15, 5);
fprintf('A 5.ª iteração chega à precisão de máquina: x_5 - 3 = %.1e (%d iterações)\n', x5 - 3, info5.nit);
Tg = [2.5000000 -3.64000 8.9120 2.9084381
      2.9084381 -0.56685 6.3898 2.9971495
      2.9971495 -0.01712 6.0114 2.9999973
      2.9999973 -0.00002 6.0000 3.0000000];
eg = abs(H(:, 6) - 3);
c = [size(H, 1) == 4 && confere(H(:, [2 6]), Tg(:, [1 4]), 7) && confere(H(:, 4), Tg(:, 2), 5) ...
       && confere(H(:, 5), Tg(:, 3), 4), ...
     confere(eg', [9.2e-2 2.9e-3 2.7e-6 2.4e-12], [3 4 7 13]), confere(H(1, 6), 2.908, 3), ...
     info.nit == 4, info.ngev == 5, info.nhev == 4, info.nfev == 0, x5 == 3];
fprintf('  tabela k = 0..3: %s | erros: %s | x1 = 2.908: %s | 4 it.: %s | n_g = 5: %s | n_H = 4: %s | n_f = 0: %s | x_5 = 3: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------- convergência quadrática
f   = @(x) x.^3 - 2*x.^2 + 4;
df  = @(x) 3*x.^2 - 4*x;
ddf = @(x) 6*x - 4;
x0 = 1;  tolg = 1e-3;
fprintf('\nConvergência quadrática: f(x) = x^3 - 2x^2 + 4, x0 = %g, tolg = %g\n', x0, tolg);
[x, ~, info] = Newton1D(f, df, ddf, x0, tolg, 20);
H = info.history;
fprintf('%3s %9s %10s %9s %9s %11s\n', 'k', 'x_k', 'df(x_k)', 'ddf(x_k)', 'x_k+1', 'df(x_k+1)');
for r = 1:size(H, 1)
  fprintf('%3d %9.5f %+10.5f %9.4f %9.5f %+11.5f\n', H(r, [1 2 4 5 6 7]));
end
e = abs([H(:, 2); x] - 4/3);
Hx = ddf(x);                         % classificar: 1 avaliação extra de f''
fprintf('|f''(x_3)| = %.1e < tolg: para em x_3 = %.5f (x* = 4/3); f''''(x_3) = %.4f > 0 (f''''(4/3) = %g): mínimo\n', ...
        abs(H(end, 7)), x, Hx, ddf(4/3));
fprintf('Erros e_k: %.4f -> %.4f -> %.4f -> %.4f;  e_3/e_2^2 = %.2f (C = 0.75)\n', e, e(4)/e(3)^2);
Tc = [1.00000 -1.00000 2.0000 1.50000  0.75000
      1.50000  0.75000 5.0000 1.35000  0.06750
      1.35000  0.06750 4.1000 1.33354  0.00081];
c = [size(H, 1) == 3 && confere(H(:, [2 4 6 7]), Tc(:, [1 2 4 5]), 5) && confere(H(:, 5), Tc(:, 3), 4), ...
     confere(abs(H(end, 7)), 8.1e-4, 5), confere(e', [0.33 0.17 0.017 0.0002], [2 2 3 4]), ...
     ddf(4/3) == 4, Hx > 0, info.ngev == 4, info.nhev == 3];
fprintf('  tabela k = 0..2: %s | |f''(x_3)| = 8.1e-4: %s | erros: %s | f''''(4/3) = 4: %s | mínimo: %s | n_g = 4: %s | n_H = 3: %s\n', ...
        simnao{c + 1});
ok = ok && all(c);
[x03, ~, i03] = Newton1D(f, df, ddf, 0.3, 1e-8, 20);
fprintf('De x0 = 0.3: x = %.6f, f''''(x) = %g < 0: máximo (ponto estacionário errado)\n', x03, ddf(x03));
c = [abs(x03) < 1e-6, ddf(x03) < 0];
fprintf('  converge para 0: %s | máximo: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------ Newton procura f' = 0 (maximizar)
f   = @(x) 8 + 5*x - 3*x.^4 - 2*x.^6;
df  = @(x) 5 - 12*(x.^3 + x.^5);
ddf = @(x) -12*(3*x.^2 + 5*x.^4);
x0 = 0.5;  tolg = 1e-7;
fprintf('\nMaximizar f(x) = 8 + 5x - 3x^4 - 2x^6, x0 = %g, tolg = %g\n', x0, tolg);
[x, fx, info] = Newton1D(f, df, ddf, x0, tolg, 20);
H = info.history;
fprintf('%3s %11s %11s %11s %11s %11s\n', 'k', 'x_k', 'f(x_k)', 'df(x_k)', 'ddf(x_k)', 'x_k+1');
for r = 1:size(H, 1)
  fprintf('%3d %11.7f %11.6f %+11.6f %11.5f %11.7f\n', H(r, 1:6));
end
Hx = ddf(x);                         % classificar: 1 avaliação extra de f''
fprintf('Após %d iterações |f''| < 1e-7: x* = %.7f, f(x*) = %.6f, f''''(x*) = %.4f < 0: máximo\n', ...
        info.nit, x, fx, Hx);
Tm = [0.5000000 10.281250  3.125000 -12.75000 0.7450980
      0.7450980 10.458621 -2.719687 -38.47906 0.6744184
      0.6744184 10.563259 -0.355311 -28.78702 0.6620756
      0.6620756 10.565490 -0.009182 -27.30912 0.6617394
      0.6617394 10.565491 -0.000007 -27.26970 0.6617391];
c = [size(H, 1) == 5 && confere(H(:, [2 6]), Tm(:, [1 5]), 7) && confere(H(:, 3:4), Tm(:, 2:3), 6) ...
       && confere(H(:, 5), Tm(:, 4), 5), ...
     info.nit == 5, confere(x, 0.6617391, 7), confere(fx, 10.565491, 6), Hx < 0];
fprintf('  tabela k = 0..4: %s | 5 it.: %s | x*: %s | f(x*): %s | máximo: %s\n', simnao{c + 1});
ok = ok && all(c);

% ------------------------------------------------ os quatro métodos
f   = @(x) x.^2 + 54./x;
df  = @(x) 2*x - 54./x.^2;
ddf = @(x) 2 + 108./x.^3;
[x5, ~, i5] = Newton1D(f, df, ddf, 1, 1e-2, 20);      % para após 5 it.
[x7, ~, i7] = Newton1D(f, df, ddf, 1, 1e-10, 20);     % para após 7 it.
fprintf('\nOs quatro métodos, Newton de x0 = 1: erro após %d it. = %.1e (n_g = %d, n_H = %d); após %d it. = %.1e\n', ...
        i5.nit, abs(x5 - 3), i5.ngev, i5.nhev, i7.nit, abs(x7 - 3));
% o slide escreve 4,5e-4
c = [i5.nit == 5, round(abs(x5 - 3)/1e-5) == 45, i5.ngev == 6, i5.nhev == 5, ...
     i7.nit == 7, floor(log10(abs(x7 - 3))) == -15];
fprintf('  5 it.: %s | erro 4,5e-4: %s | n_g = 6: %s | n_H = 5: %s | 7 it.: %s | erro da ordem de 1e-15: %s\n', simnao{c + 1});
ok = ok && all(c);

fprintf('\nconfere com os slides: %s\n', simnao{ok + 1});
